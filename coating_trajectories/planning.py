from numpy import sqrt, dot, concatenate, arange, array
from numpy import transpose, linalg, sum, cross, zeros, eye
from openravepy import IkFilterOptions, Ray
from math import pi, cos, sin, atan2
from scipy.optimize import minimize, linprog
from copy import copy
from path_filters import filter_trajectories
from mathtools import central_difference
import mathtools
import time
import logging
from openravepy import databases, IkParameterization
from profilestats import profile

"""
Main package for robot joints' positions and velocities planning,
robot base calculation, torque and manipulability analysis.
"""

def compute_angular_velocities(turbine, joints_trajectory, trajectory_index):

    if (trajectory_index>2) and ((len(joints_trajectory)-trajectory_index)>3):
        return central_difference(turbine, joints_trajectory, trajectory_index)
    else: return None

def torque_computation(turbine, joints, w, alpha):
    # INVERSE DYNAMICS ComputeInverseDynamics
    with turbine.robot:
        turbine.robot.SetDOFValues(joints, turbine.robot.GetActiveDOFIndices(), turbine.robot.CheckLimitsAction.CheckLimitsThrow)
        turbine.robot.SetDOFVelocities(w, turbine.robot.GetActiveDOFIndices(), turbine.robot.CheckLimitsAction.CheckLimitsThrow)
        torques = turbine.robot.ComputeInverseDynamics(alpha)

        if any(torques < turb.robot.GetDOFMaxTorque):
            return False#THROW EXCEPTION
        else:
            return True
    return torques

def sensibility(turbine, point_normal, w, alpha):
    
    normal_vec = point_normal[3:6]
    tangent_vec = cross(point_normal[0:3],normal_vec)
    tangent_vec = tangent_vec/sqrt(dot(tangent_vec,tangent_vec))
    perp_vec = cross(normal_vec,tangent_vec)
    
    Jpos = turbine.manipulator.CalculateJacobian()
    Hpos = turbine.robot.ComputeHessianTranslation(turbine.manipulator.GetArmDOF(),
                                                   turbine.manipulator.GetEndEffectorTransform()[0:3,3])
    Hpos = dot(Hpos,w)
    theta_limits = zip(-2*turbine.robot.GetDOFResolutions(),2*turbine.robot.GetDOFResolutions())
    w_limits = zip(-w*0.01,w*0.01) #HARDCODED 1% error
    limts = tuple(theta_limits + w_limits)

    Hpos_tan = dot(Hpos,tangent_vec)
    Jpos_tan = dot(tangent_vec, Jpos)
    errorgain_tan = contcatenate((Jpos_tan,Hpos_tan))
    velocity_tan_error = (linprog(errorgain_tan, bounds=limts).get('fun'),-linprog(-errorgain_tan, bounds=limts).get('fun'))

    Jpos_normal = dot(normal_vec, Jpos)
    position_normal_error = (linprog(Jpos_normal, bounds=theta_limits).get('fun'),-linprog(-Jpos_normal, bounds=theta_limits).get('fun'))

    Jpos_perp = dot(perp_vec, Jpos)
    position_perp_error = (linprog(Jpos_perp, bounds=theta_limits).get('fun'),-linprog(-Jpos_perp, bounds=theta_limits).get('fun'))

    Jw = turbine.manipulator.CalculateAngularVelocityJacobian()
    x_dir = turbine.manipulator.GetTransform()[0:3,0]
    nhat = mathtools.hat(normal_vec)
    xn = dot(x_dir,nhat)
    Jcos = -dot(xn,Jw)
    cos = -dot(x_dir,normal_vec)
    dcos = (linprog(Jcos, bounds=theta_limits).get('fun'),-linprog(-Jcos, bounds=theta_limits).get('fun'))
    angle_error = tuple(np.arccos(cos+dcos))

    return velocity_tan_error, position_normal_error, position_perp_error, angle_error
    
def compute_robot_joints(turbine, trajectory, trajectory_index):
    """
    Iterates points of the trajectories, computing optimal robot's joints
    (minimizing orientation error).

    Keyword arguments:
    turbine -- turbine object
    trajectory -- trajectory to coat
    trajectory_index -- where to begin. Index of a feasible point in the trajectory.
    joint_solutions -- initial joint_solutions
    """

    joint_solutions = []
    robot = turbine.robot
    logging.info('Computing robot joints from point: ' + str(trajectory[trajectory_index]))
    
    # Compute inverse kinematics and find a solution
    ans, joint_solution, _ = compute_feasible_point(turbine, trajectory[trajectory_index])
    
    if not ans:
        return []

    joint_solutions.append(joint_solution)
    # Find solutions for next points
    robot.SetDOFValues(joint_solution)
    for index in range(trajectory_index+1, len(trajectory)):
        res = orientation_error_optimization(turbine, trajectory[index])
        if trajectory_constraints(turbine, res, trajectory[index]):
            joint_solutions.append(res.x)
            robot.SetDOFValues(res.x)
        else:
            ans, joint_solution, _ = compute_feasible_point(turbine, trajectory[index])
            if not ans:
                logging.info('--NOT ANS--')
                return joint_solutions
            robot.SetDOFValues(joint_solution)
            previous_joint_solutions = backtrack(turbine, trajectory[:index+1])
            if len(previous_joint_solutions)==0:
                logging.info('----previous_joint_solutions----')
                return joint_solutions
            joint_solutions = previous_joint_solutions
            robot.SetDOFValues(joint_solution)
    return joint_solutions

#@profile(print_stats=10, dump_stats=True)
def compute_first_feasible_point(turbine, trajectory):
    """
    Method to compute the first feasible point in the trajectory: where to start.

    Keyword arguments:
    turbine -- turbine object
    trajectory -- trajectory to coat
    """

    robot = turbine.robot
    with robot:
        for i in range(0,len(trajectory)):
            ans, sol, tolerance = compute_feasible_point(turbine, trajectory[i])
            if ans:
                logging.info('Feasible point found: '+str(trajectory[i]))
                return i, sol, tolerance
        logging.info('Feasible point not found.')
        raise ValueError('No solution for given trajectory')

def compute_feasible_point(turbine, point):
    """
    Method to compute if point is feasible.

    Keyword arguments:
    turbine -- turbine object
    point -- point to coat
    """

    with turbine.robot:
        iksol, tolerance = ik_angle_tolerance(turbine, point)
        if len(iksol)>0:
            logging.info('Point with tolerance found: '+str(point))
            turbine.robot.SetDOFValues(best_joint_solution_regarding_manipulability(iksol, turbine.robot))
            res = orientation_error_optimization(turbine, point)
            if(trajectory_constraints(turbine, res, point)):
                return True, res.x, tolerance
        return False, [], []

def orientation_error_optimization(turbine, point, tol=1e-3):
    """
    Minimizes the orientation error, given an initial configuration of the robot
    and the desired point to coat (with its normal vector). The constraints of the minimi-
    zation are: to reach de point position (equality).

    Keyword arguments:
    robot -- the robot.
    point -- 6x1 vector is the goal (point to be coated). (x,y,z) cartesian position
    and (nx,ny,nz) is the normal vector of the point, w.r.t. the world frame.
    """
    
    robot = turbine.robot
    q0 = robot.GetDOFValues()
    manip = robot.GetActiveManipulator()
    lower_limits, upper_limits = robot.GetActiveDOFLimits()
    lower_limits = lower_limits+0.01
    upper_limits = upper_limits-0.01
    
    with robot:
        def func(q):
            robot.SetDOFValues(q, manip.GetArmIndices())
            Rx = manip.GetTransform()[0:3,0]
            Rx = Rx/linalg.norm(Rx)
            return dot(point[3:6], Rx)
        
        def position_cons(q):
            robot.SetDOFValues(q)
            pos = manip.GetTransform()[0:3,3]
            dif = pos-point[0:3]
            return dot(dif,dif)/tol
        
        cons = ({'type':'eq', 'fun': position_cons})
        
        bnds = tuple([(lower_limits[i],upper_limits[i]) for i in range(0,len(lower_limits))])
        
        res = minimize(func, q0, constraints=cons, method='SLSQP',
                       bounds = bnds, options={'disp': False})
    return res    

def orientation_cons(turbine, point):
    """
    Return True if the orientation error is inside limits.
    """

    Rx = turbine.robot.GetActiveManipulator().GetTransform()[0:3,0]
    Rx = Rx/linalg.norm(Rx)
    return (dot(point[3:6], Rx)) + cos(turbine.config.coating.angle_tolerance) <= 0

def backtrack(turbine, trajectory):
    """
    Call when optimization fails. Previous solution should be recomputed.

    Keyword arguments:
    turbine -- turbine object
    trajectory -- already computed coated points.
    """

    logging.info('Starting backtracking at point: ' + str(trajectory[-1]))
    robot = turbine.robot
    Q = []
    with robot:
        for point in reversed(trajectory):
            res = orientation_error_optimization(turbine, point, tol=1e-3)
            if trajectory_constraints(turbine, res, point):
                Q.append(res.x)
                robot.SetDOFValues(res.x)
            else: return []
    return list(reversed(Q))


def ik_angle_tolerance(turbine, point):
    """
    Solve the inverse kinematics given point (IKFast) with maximum tolerance angle.

    Keyword arguments:
    turbine -- turbine object
    point -- point to be coated is a 6D array, which (x,y,z) cartesian position
    and (nx,ny,nz) is the normal vector of the point, w.r.t. the world frame.
    """

    robot = turbine.robot
    manip = robot.GetActiveManipulator()
   
    iksol = ikfast(robot, point)
    if len(iksol)>0:
        return iksol, []
    
    # Compute solution with maximum distance and angle tolerances
    angle = turbine.config.coating.angle_tolerance
    number_of_angles = 72 #HARDCODED 

    normal_tol = dot(point[3:6],transpose(mathtools.Raxis(
        mathtools.compute_perpendicular_vector(point[3:6]), angle)))
    normal_tol = normal_tol/linalg.norm(normal_tol)
    
    k = 2*pi/number_of_angles
    counter = 0
    for i in range(0, number_of_angles):
        alfa = k*i
        iksol = ikfast(robot, concatenate((point[0:3],dot(
            normal_tol, transpose(mathtools.Raxis(point[3:6],alfa))))))
        if len(iksol)>0:
            counter+=1
            if counter==2:
                return iksol, ['angle', angle]
    return [], []

def ik_angle_tolerance_normal_plane(turbine, point, iter_surface, angle_discretization = 5*pi/180):
    """
    Solve the inverse kinematics given point (IKFast) with maximum tolerance angle
    discretizing only on the normal plane.
    Return example:
    iksol = [joint_solution_1, joint_solution_2, ..., joint_solution_n]
    angle_tolerance = [30*pi/180, 30*pi/180, ..., 0*pi/180, ..., -30*pi/180]
    joint_solution_1.shape = (1,6)

    Keyword arguments:
    turbine -- turbine object
    point -- point to be coated is a 6D array, which (x,y,z) cartesian position
    and (nx,ny,nz) is the normal vector of the point, w.r.t. the world frame.
    iter_surface -- surface used to generate the trajectories to compute the point
    tangent.
    angle_discretization -- in radians.
    """

    iksol = []
    angle_tolerance = []
    
    robot = turbine.robot
    
    if angle_discretization == 0:
        sols = ikfast(robot, point)
        for sol in sols:
            iksol.append(sol)
            angle_tolerance.append(angle)
        return iksol, angle_tolerance
    
    manip = robot.GetActiveManipulator()
    tolerance = turbine.config.coating.angle_tolerance
   
    tan = mathtools.surfaces_tangent(point, iter_surface)

    angle_discretization = min(abs(angle_discretization),
                               tolerance)
    
    for angle in arange(-tolerance, tolerance + angle_discretization, angle_discretization):
        normal_tol = dot(point[3:6],transpose(mathtools.Raxis(tan, angle)))
        sols = ikfast(robot, concatenate((point[0:3], normal_tol)))
        for sol in sols:
            iksol.append(sol)
            angle_tolerance.append(angle)
            
    return iksol, angle_tolerance 

def Fcost(turbine, start, actual_ik_cells, actual_tolerance_cells, goal):

    robot = turbine.robot

    actual_ik_cells = array(actual_ik_cells)
    actual_tolerance_cells = array(actual_tolerance_cells)   
    delta_t = turbine.config.model.trajectory_step/turbine.config.coating.coating_speed

    # Verifying joint velocities consistency from start and goal
    joint_distance_from_start = actual_ik_cells - start
    joint_distance_from_goal = actual_ik_cells - goal
    ws_in_limits_start = [w.all() for w in ((abs(joint_distance_from_start/delta_t) - robot.GetDOFVelocityLimits()) < 0)]
    ws_in_limits_goal = [w.all() for w in ((abs(joint_distance_from_goal/delta_t) - robot.GetDOFVelocityLimits()) < 0)]
    ws_in_limits = ws_in_limits_start and ws_in_limits_goal
    
    actual_ik_cells = actual_ik_cells[ws_in_limits]
    actual_tolerance_cells = actual_tolerance_cells[ws_in_limits]
    joint_distance_from_start = joint_distance_from_start[ws_in_limits]
    joint_distance_from_goal = joint_distance_from_goal[ws_in_limits]

    # Computing distance from start and goal
    if (len(joint_distance_from_goal) > 0) and (len(joint_distance_from_start) > 0):
        G_cost = sqrt(sum(joint_distance_from_start*joint_distance_from_start,1))
        H_cost = sqrt(sum(joint_distance_from_goal*joint_distance_from_goal,1))
        # Heuristic cost
        F_cost = G_cost + H_cost + actual_tolerance_cells
        return F_cost, actual_ik_cells
    else:
        return [], []

def Gcost(turbine, start, actual_ik_cells, actual_tolerance_cells):

    robot = turbine.robot

    actual_ik_cells = array(actual_ik_cells)
    actual_tolerance_cells = array(actual_tolerance_cells)   
    delta_t = turbine.config.model.trajectory_step/turbine.config.coating.coating_speed

    # Verifying joint velocities consistency from start and goal
    joint_distance_from_start = actual_ik_cells - start
    ws_in_limits = [w.all() for w in ((abs(joint_distance_from_start/delta_t) - robot.GetDOFVelocityLimits()) < 0)]
    
    actual_ik_cells = actual_ik_cells[ws_in_limits]
    actual_tolerance_cells = actual_tolerance_cells[ws_in_limits]
    joint_distance_from_start = joint_distance_from_start[ws_in_limits]

    # Computing distance from start
    if len(joint_distance_from_start) > 0:
        G_cost = sum(abs(joint_distance_from_start/robot.GetDOFVelocityLimits()),1)
        return G_cost, actual_ik_cells, actual_tolerance_cells
    else:
        return [], [], []

def ikfast(robot, point):
    """
    Call openrave IKFast. It computes the inverse kinematic for the point.
    It returns all solutions.

    keyword arguments:
    robot -- the robot. 
    point -- point to coat is a 6D array, which (x,y,z) cartesian position
    and (nx,ny,nz) is the normal vector of the point, w.r.t. the world frame.
    """
    
    with robot:
        manip = robot.GetActiveManipulator()
        Tee = manip.GetTransform()
        Rx = Tee[0:3,0]
        Rx = Rx/linalg.norm(Rx)
        Rab = mathtools.Rab(Rx, -point[3:6])
        
        T = eye(4)
        T[0:3,0:3] = dot(Rab,Tee[0:3,0:3])
        T[0:3,3] = point[0:3]
        solutions = robot.GetActiveManipulator().FindIKSolutions(T, IkFilterOptions.CheckEnvCollisions)

        if len(solutions)>0:
            if len(solutions.shape)==1:
                solutions = solutions.reshape((1,solutions.shape[0]))
                
        return solutions

def ikfast_position(robot, point):
    """
    Call openrave IKFast Translation3D. It computes the inverse kinematic for the point.
    It returns one solution.

    keyword arguments:
    robot -- the robot. 
    point -- point to coat is a 6D array, which (x,y,z) cartesian position
    and (nx,ny,nz) is the normal vector of the point, w.r.t. the world frame.
    """

    with robot:
        ikparam = IkParameterization(point[0:3],IkParameterization.Type.Translation3D)
        sol = robot.GetActiveManipulator().FindIKSolution(ikparam, IkFilterOptions.CheckEnvCollisions)
        return sol, []

def ikfast_5d(robot, point):
    """
    Call openrave IKFast TranslationDirection5D. It computes the inverse kinematic for the point.
    It returns all solutions.

    keyword arguments:
    robot -- the robot. 
    point -- point to coat is a 6D array, which (x,y,z) cartesian position
    and (nx,ny,nz) is the normal vector of the point, w.r.t. the world frame.
    """

    with robot:
        ikparam = IkParameterization(Ray(point[0:3],-point[3:6]),IkParameterization.Type.TranslationDirection5D)
        solutions = robot.GetActiveManipulator().FindIKSolutions(ikparam, IkFilterOptions.CheckEnvCollisions)
        return sol, []

def check_dof_limits(robot, q):
    """
    Method checks if the joints limits are inside limits with a safe tolerance.
    return True if limits were respected.
    keyword arguments:
    robot -- the robot. 
    q -- DOF values.
    """
    joints = robot.GetJoints()
    for i in range(0,len(q)):
        l,u = joints[i].GetLimits()
        l = l[0]
        u = u[0]
        if not q[i]>=l+0.01 or not q[i]<=u-0.01:
            return False
    return True

def compute_manipulability_det(robot, joint_configuration):
    """
    Compute robot manipulability as described in:
    Yoshikawa, 'Manipulability of robotic mechanisms',
    International Journal of Robotics Research, vol. 4
    no. 2, pp. 3-9, 1985.
    
    keyword arguments:
    robot -- the robot. 
    joint_configuration -- q Nx1, joint configuration.
    """
                                 
    manip = robot.GetActiveManipulator()
    with robot:
        robot.SetDOFValues(joint_configuration)
        Jpos = manip.CalculateJacobian()
        Jori = manip.CalculateAngularVelocityJacobian()
        J = concatenate((Jpos,Jori))
    return sqrt(linalg.det(dot(transpose(J),J))), sqrt(linalg.det(dot(Jpos,transpose(Jpos)))), sqrt(linalg.det(dot(Jori,transpose(Jori))))

def best_joint_solution_regarding_manipulability(joint_solutions, robot):
    """
    Given a list of joint solutions for a specific point, the function computes
    maniulability and returns the best solution regarding manipulability criteria.
    """
                                 
    biggest_manipulability = 0
    temp_q = []
    for q in joint_solutions:
        try:
            manipulability, _, _ = compute_manipulability_det(robot, q)
            if manipulability>biggest_manipulability:
                biggest_manipulability = manipulability
                temp_q = q
        except IndexError:
            raise IndexError('There is no solution for the given trajectory.')
    return temp_q

def trajectory_constraints(turbine, res, point):
    """
    Check robot self and environment collision, optimization result,
    and angle tolerance.
    """

    robot = turbine.robot

    logging.basicConfig(filename='trajectory_constraints.log', level=logging.DEBUG)
    # Verifying optimization solution
    if not res.success:
        logging.info('Optimization failed in point: '+str(point)+', due to optimization fail.')
        return False

    # Verifying environment and self collision
    with robot:
        robot.SetDOFValues(res.x)
        if turbine.check_robot_collision():
            logging.info('Trajectory terminates in point: '+str(point)+', due to env collision detection.')
            return False
        if robot.CheckSelfCollision():
            logging.info('Trajectory terminates in point: '+str(point)+', due to self collision detection.')
            return False

        # Verifying angle tolerance
        if not orientation_cons(turbine, point):
            logging.info('Trajectory terminates in point: '+str(point)+', due to orientation constraint.')
            return False
    return True

def set_ikmodel_translation3d(robot):
    ikmodel = databases.inversekinematics.InverseKinematicsModel(
        robot=robot, iktype=IkParameterization.Type.Translation3D)
    if not ikmodel.load():
        ikmodel.autogenerate()
    return

def set_ikmodel_transform6D(robot):
    ikmodel = databases.inversekinematics.InverseKinematicsModel(
        robot=robot, iktype=IkParameterization.Type.Transform6D)
    if not ikmodel.load():
        ikmodel.autogenerate()
    return

def set_ikmodel_translationdirection5d(robot):
    ikmodel = databases.inversekinematics.InverseKinematicsModel(
        robot=robot, iktype=IkParameterization.Type.TranslationDirection5D)
    if not ikmodel.load():
        ikmodel.autogenerate()
    return
