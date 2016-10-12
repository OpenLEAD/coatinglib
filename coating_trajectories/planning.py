from numpy import sqrt, dot, concatenate, arange, array, zeros, transpose, linalg, sum, cross
from openravepy import IkFilterOptions
from math import pi, cos, sin, atan2
from scipy.optimize import minimize
from copy import copy
from collections import deque
import mathtools
import time
import logging

"""
Main package for robot joints' positions and velocities planning,
robot base calculation, torque and manipulability analysis.
"""

def _std_robot_filter(turbine, trajectories):
    raise ValueError("No trajectory filter for "+turbine.robot.GetName()+" robot. Create new function.")

def mh12_filter(turbine, trajectories):
    pistol = 0.3
    flame = 0.23
    _working_radius_squared = (1.285+pistol+flame)**2
    def distance_robot_squared(trajectory, robot):
        delta = robot.GetJoints()[1].GetAnchor() - trajectory
        return sum(delta*delta,1)

    filtered_trajectories = []
    for trajectory in trajectories:
        trajectory = array(trajectory)
        filtered_trajectory = trajectory[distance_robot_squared(trajectory[:,0:3], turbine.robot) < _working_radius_squared]
        if len(filtered_trajectory)>0:
            filtered_trajectories.append(filtered_trajectory)

    return filtered_trajectories


_filter_options = {'mh12': mh12_filter}    

def filter_trajectories(turbine, trajectories):
    name = turbine.robot.GetName()
    _filter_options.get(name,_std_robot_filter)(turbine, trajectories)

def workspace_limit(turbine,point_normal, joints_trajectory, trajectory_index):

    # 'j' for joints, so it doesnt clumpsy the equations
    j = deque(joints_trajectory)
    j.rotate(-trajectory_index)

    h = turbine.config.model.trajectory_step
    
    # Joints velocity - Central Difference (h**6 order error)
    w = ( (j[3]-j[-3]) + 9*(j[-2]-j[2]) + 45*(j[1]-j[-1]) )/(60.0*h)

    # Joints acceleration - Central Difference (h**6 order error)
    alpha = ( 2*(j[-3]+j[3]) - 27*(j[-2]+j[2]) + 270*(j[-1]+j[1]) - 490*j[0] )/(180.0*h**2)

    # INVERSE DYNAMICS ComputeInverseDynamics
    with turbine.robot:
        turbine.robot.SetDOFValues(j[0], turb.robot.GetActiveDOFIndices(), turbine.robot.CheckLimitsAction.CheckLimitsThrow)
        turbine.robot.SetDOFVelocities(w, turb.robot.GetActiveDOFIndices(), turbine.robot.CheckLimitsAction.CheckLimitsThrow)
        torques = turb.robot.ComputeInverseDynamics(alpha)
        if any(torques < turb.robot.GetDOFMaxTorque):
            return False#THROW EXCEPTION
        else:
            return True

    #Frame_transform=np.transpose(turb.manipulator.GetTransform()[0:3,0:3])
    normal_vec = point_normal[3:6]
    #normal_vec = np.dot(R,normal_vec)
    tangent_vec = cross(point_normal[0:3],normal_vec)
    tangent_vec = tangent_vec/sqrt(dot(tangent_vec,tangent_vec))
    #tangent_vec = np.dot(R,tangent_vec)
    perp_vec = cross(normal_vec,tangent_vec)
    #perp_vec = np.dot(R,perp_vec)

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
    position_normal_error = (linprog(Jpos_normal, bounds=limts).get('fun'),-linprog(-Jpos_normal, bounds=limts).get('fun'))

    Jpos_perp = dot(perp_vec, Jpos)
    position_perp_error = (linprog(Jpos_perp, bounds=limts).get('fun'),-linprog(-Jpos_perp, bounds=limts).get('fun'))

    return w,alpha,torques,velocity_tan_error,position_normal_error,position_perp_error 
    

def compute_robot_joints(turbine, trajectory, trajectory_index, joint_solutions = []):
    """
    Iterates points of the trajectories, computing optimal robot's joints
    (minimizing orientation error).

    Keyword arguments:
    turbine -- turbine object
    trajectory -- trajectory to coat
    trajectory_index -- where to begin
    joint_solutions -- initial joint_solutions
    """

    robot = turbine.robot

    # Compute inverse kinematics and find a solution
    joint_solution, _ = inverse_kinematics(turbine, trajectory[trajectory_index])   
    best_joint_solution = best_joint_solution_regarding_manipulability(joint_solution, robot)
    robot.SetDOFValues(best_joint_solution)

    # Find solutions for previous points
    previous_joint_solutions = backtrack(turbine, trajectory[:trajectory_index+1])
    if len(previous_joint_solutions)==0:
        return joint_solution

    if len(joint_solutions)>0: joint_solutions = previous_joint_solutions + joint_solutions
    else: joint_solutions = previous_joint_solutions

    # Find solutions for next points
    for index in range(trajectory_index+1, len(trajectory)):
        res = orientation_error_optimization(turbine, trajectory[index])
        if trajectory_constraints(turbine, res, trajectory[index]):
            joint_solutions.append(res.x)
        else:
            return compute_robot_joints(turbine, trajectory, index, joint_solutions)
    return joint_solutions

def orientation_error_optimization(turbine, point, tol=1e-3):
    """ Minimizes the orientation error, given an initial configuration of the robot
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
            T = manip.GetTransform()
            Rx = T[0:3,0]
            Rx = Rx/sqrt(dot(Rx,Rx))
            return dot(point[3:6], Rx)
        
        def position_cons(q):
            robot.SetDOFValues(q, manip.GetArmIndices())
            T = manip.GetTransform()
            pos = T[0:3,3]
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
    
    cos_tolerance = cos(turbine.config.coating.angle_tolerance)
    robot = turbine.robot
    manip = robot.GetActiveManipulator()
    T = manip.GetTransform()
    Rx = T[0:3,0]
    Rx = Rx/sqrt(dot(Rx,Rx))
    return (-dot(point[3:6], Rx) - cos_tolerance)>=0

def backtrack(turbine, trajectory):
    """
    Call when optimization fails. Previous solution should be recomputed.

    Keyword arguments:
    turbine -- turbine object
    trajectory -- already computed coated points.
    """

    robot = turbine.robot
    Q = []
    
    for point in reversed(trajectory):
        res = orientation_error_optimization(turbine, point, tol=1e-3)
        if trajectory_constraints(turbine, res, point):
            Q.append(res.x)
        else: return []
    return list(reversed(Q))

def inverse_kinematics(turbine, point):
    """
    Solve the inverse kinematics given point (IKFast).

    Keyword arguments:
    turbine -- turbine object
    point -- point to be coated is a 6D array, which (x,y,z) cartesian position
    and (nx,ny,nz) is the normal vector of the point, w.r.t. the world frame.
    """

    robot = turbine.robot
    manip = robot.GetActiveManipulator()
    coating_tolerance_max = turbine.config.coating.max_distance
    coating_tolerance_min = turbine.config.coating.min_distance

    # Compute solution without tolerances
    iksol = ikfast(robot, point)
    if len(iksol)>0:
        return iksol, []

    # Compute solution with distance and angle tolerances
    angles = []
    if turbine.config.coating.angle_tolerance>0:
        angles = arange(0, turbine.config.coating.angle_tolerance, 0.001)
    number_of_angles = 10

    distance_tolerance = arange(coating_tolerance_min, coating_tolerance_max, 1e-3)
    if len(distance_tolerance)>0:
        for distance in arange(coating_tolerance_min, coating_tolerance_max, 1e-3):
            new_point = concatenate((point[0:3]+distance*point[3:6], point[3:6]))
            iksol = ikfast(robot, new_point)
            if len(iksol)>0: return iksol, ['distance', distance] 
            if len(angles)>0:
                for angle in angles:
                    Rz = mathtools.Raxis(mathtools.compute_perpendicular_vector(manip.GetDirection()),
                                         angle)
                    p = dot(manip.GetDirection(), transpose(Rz))
                    k = 2*pi/number_of_angles
                    for i in range(0, number_of_angles):
                        alfa = k*i
                        Rd = mathtools.Raxis(manip.GetDirection(),alfa)
                        iksol = ikfast(robot,
                                       concatenate((point[0:3],dot(p, transpose(Rd)))))
                        if len(iksol)>0:return iksol, ['angle', angle]

    # No distance tolerance, but angle tolerance
    else:
        if len(angles)>0:
            for angle in angles:
                Rz = mathtools.Raxis(mathtools.compute_perpendicular_vector(manip.GetDirection()),
                                     angle)
                p = dot(manip.GetDirection(), transpose(Rz))
                k = 2*pi/number_of_angles
                for i in range(0, number_of_angles):
                    alfa = k*i
                    Rd = mathtools.Raxis(manip.GetDirection(),alfa)
                    iksol = ikfast(robot,
                                   concatenate((point[0:3],dot(p, transpose(Rd)))))
                    if len(iksol)>0:return iksol, ['angle', angle]
    return [], []
            
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
        T = zeros((4,4))
        T[0:3,0:3] = mathtools.Rab(manip.GetDirection(), point[3:6])
        T[0:3,3] = point[0:3]
        T[3,0:4] = [0,0,0,1]
        solutions = robot.GetActiveManipulator().FindIKSolutions(T, IkFilterOptions.CheckEnvCollisions)

        if len(solutions)>0:
            if len(solutions.shape)==1:
                solutions = solutions.reshape((1,solutions.shape[0]))
                
        return solutions 

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
    Given a list of joint solutions for a specific point, the function computes maniulability and returns the
    best solution regarding manipulability criteria.
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
        logging.info('Optimization failed in point: '+str(index)+'/'+str(len(trajectory))+', due to optimization fail.')
        return False

    # Verifying environment and self collision
    robot.SetDOFValues(res.x)
    if turbine.check_robot_collision():
        logging.info('Trajectory terminates in point: '+str(index)+'/'+str(len(trajectory))+', due to env collision detection.')
        return False
    if robot.CheckSelfCollision():
        logging.info('Trajectory terminates in point: '+str(index)+'/'+str(len(trajectory))+', due to self collision detection.')
        return False

    # Verifying angle tolerance
    if not orientation_cons(turbine, point):
        logging.info('Trajectory terminates in point: '+str(index)+'/'+str(len(trajectory))+', due to orientation constraint.')
        return False

    return True
