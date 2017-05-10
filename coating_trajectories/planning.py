from numpy import sqrt, dot, concatenate, arange, array
from numpy import abs, zeros, cumsum, minimum
from numpy import transpose, linalg, sum, cross, zeros, eye, max, inf
from numpy import arccos, maximum, random
from numpy.linalg import norm
from openravepy import IkFilterOptions, interfaces, databases
from openravepy import IkParameterization
from openravepy import CollisionOptions, RaveCreateCollisionChecker, CollisionReport
from math import pi, cos, sin, atan2
from scipy.optimize import minimize, linprog
from copy import copy
from path_filters import filter_trajectories
from mathtools import central_difference
import mathtools
import dijkstra2

"""
Main package for robot joints' positions and velocities planning,
robot base calculation, torque and manipulability analysis.
"""

def compute_angular_velocities(turbine, joints_trajectory, trajectory_index):

    if (trajectory_index>2) and ((len(joints_trajectory)-trajectory_index)>3):
        return central_difference(turbine, joints_trajectory, trajectory_index)
    else: return None
    
def compute_general_velocities(turbine, joints_trajectory, points):
    """
    points = rays[:,0:3]
    0)Compute times based on h and distances calculated from rays
    
    Loop joints_trajectory/rays
    1)Take at most 7 points of the rays per joints, the optimal case is centering the actual ray
    2)Order elements for best approximation like: 0,1,-1,2,-2..
    3)Use general_finite_difference(time, joints, times) to compute difference, it outputs the estimated joint, w , alpha
    """
    h = turbine.config.coating.coating_speed
    dpoints = array(points[1:])-array(points[:-1])
    times = cumsum(array([0.]+list(norm(dpoints,axis=1)))/h)

    N = len(times)

    w_list = []
    alpha_list = []
    
    for i,time in enumerate(times):
        past = range(i)[::-1]
        future = range(i,N)
        order = list(mathtools.roundrobin(future,past))
        order = order[:min(len(order),7)]
        _,w,alpha = mathtools.general_finite_difference(time, array(joints_trajectory)[order], times[order])
        w_list += [w]
        alpha_list += [alpha]
    return w_list, alpha_list, times

def torque_computation(turbine, joints, w, alpha, verify = False):
    """
    Apply joints position 'joints', joints velocities 'w' and joints acceleration 'alpha'
    Return torques - (6,) array
    verify - if True function returns False when a torque is bigger than maximum
    """
    with turbine.robot:
        gun = turbine.robot.GetLink('Gun')
        flame_force = {gun.GetIndex(): list(-gun.GetTransform()[0:3,2]*turbine.config.coating.flame_thrust)+[0,0,0]} 
        turbine.robot.SetDOFValues(joints, turbine.robot.GetActiveDOFIndices(), turbine.robot.CheckLimitsAction.CheckLimitsSilent)
        turbine.robot.SetDOFVelocities(w,turbine.robot.CheckLimitsAction.CheckLimitsSilent,turbine.robot.GetActiveDOFIndices())
        torques = turbine.robot.ComputeInverseDynamics(alpha,flame_force)

    if any(torques > turbine.robot.GetDOFMaxTorque) and verify:
        return False
    else:
        return torques

def sensibility(turbine, ray, w, alpha):
    
    normal_vec = ray[3:6]
    tangent_vec = cross(ray[0:3],normal_vec)
    tangent_vec = tangent_vec/sqrt(dot(tangent_vec,tangent_vec))
    perp_vec = cross(normal_vec,tangent_vec)
    
    Jpos = turbine.manipulator.CalculateJacobian()
    Hpos = turbine.robot.ComputeHessianTranslation(turbine.manipulator.GetArmDOF(),
                                                   turbine.manipulator.GetEndEffectorTransform()[0:3,3])
    Hpos = dot(Hpos,w)
    theta_limits = zip(-turbine.robot.GetDOFResolutions(),turbine.robot.GetDOFResolutions())
    w_limits = zip(-abs(w)*0.001,abs(w)*0.001) #HARDCODED 1% error
    limits = tuple(theta_limits + w_limits)

    Hpos_tan = dot(Hpos,tangent_vec)
    Jpos_tan = dot(tangent_vec, Jpos)
    errorgain_tan = concatenate((Jpos_tan,Hpos_tan))
    velocity_tan_error = (linprog(errorgain_tan, bounds=limits).get('fun'),-linprog(-errorgain_tan, bounds=limits).get('fun'))

    Jpos_normal = dot(normal_vec, Jpos)
    position_normal_error = (linprog(Jpos_normal, bounds=theta_limits).get('fun'),-linprog(-Jpos_normal, bounds=theta_limits).get('fun'))

    Jpos_perp = dot(perp_vec, Jpos)
    position_perp_error = (linprog(Jpos_perp, bounds=theta_limits).get('fun'),-linprog(-Jpos_perp, bounds=theta_limits).get('fun'))

    Jw = turbine.manipulator.CalculateAngularVelocityJacobian()
    x_dir = turbine.manipulator.GetTransform()[0:3,0]
    nhat = mathtools.hat(normal_vec)
    xn = dot(x_dir,nhat)
    Jcos = -dot(xn,Jw)
    cosn = -dot(x_dir,normal_vec)
    dcosn = (linprog(Jcos, bounds=theta_limits).get('fun'),-linprog(-Jcos, bounds=theta_limits).get('fun'))
    angle_error = tuple(arccos(minimum(cosn+dcosn,1.)) - arccos(cosn))

    return velocity_tan_error, position_normal_error, position_perp_error, angle_error
    
def compute_robot_joints(turbine, trajectory, trajectory_index, iter_surface = None):
    """
    Iterates points of the trajectories, computing optimal robot's joints
    (minimizing orientation error). For the first point, it uses ikfast.

    Keyword arguments:
    turbine -- turbine object
    trajectory -- trajectory to coat
    trajectory_index -- where to begin. Index of a feasible point in the trajectory.
    joint_solutions -- initial joint_solutions
    """

    joint_solutions = []
    robot = turbine.robot
    
    # Compute inverse kinematics and find a solution
    iksol, tolerance = compute_feasible_point(turbine, trajectory[trajectory_index], iter_surface)
    if len(iksol)>0:
        sol = best_joint_solution_regarding_manipulability(
            iksol, tolerance, turbine.robot)
        robot.SetDOFValues(sol)
        if(trajectory_constraints(turbine, sol, trajectory[trajectory_index])):
            joint_solutions.append(sol)
        else: return []
    else: return []

    # Find solutions for next points
    for index in range(trajectory_index+1, len(trajectory)):
        res = orientation_error_optimization(turbine, trajectory[index])
        if trajectory_constraints(turbine, res.x, trajectory[index]):
            joint_solutions.append(res.x)
            robot.SetDOFValues(res.x)
        else:
            return joint_solutions
    return joint_solutions

def compute_robot_joints_opt(turbine, trajectory, trajectory_index, iter_surface = None):
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
    iksol, _ = compute_feasible_point(turbine, trajectory[trajectory_index], iter_surface)
    best_joint_solutions = []

    for ik in iksol:
        robot.SetDOFValues(ik)
        joint_solutions = []
        joint_solutions.append(ik)
        for index in range(trajectory_index+1, len(trajectory)):
            res = orientation_error_optimization(turbine, trajectory[index])
            if res.success:
                if trajectory_constraints(turbine, res.x, trajectory[index]):
                    joint_solutions.append(res.x)
                    robot.SetDOFValues(res.x)
                else: break
            else: break
        else:
            return joint_solutions
        if len(joint_solutions)>len(best_joint_solutions):
            best_joint_solutions = joint_solutions
    return best_joint_solutions


def compute_first_feasible_point(turbine, trajectory, iter_surface = None):
    """
    Method to compute the first feasible point in the trajectory: where to start.

    Keyword arguments:
    turbine -- turbine object
    trajectory -- trajectory to coat
    """

    robot = turbine.robot
    with robot:
        for i in range(0,len(trajectory)):
            sols, tolerances = compute_feasible_point(turbine, trajectory[i], iter_surface)
            if len(sols)>0:
               return i, sols, tolerances
        raise ValueError('No solution for given trajectory')

def compute_feasible_point(turbine, point, iter_surface = None):
    """
    Method to compute if point is feasible.

    Keyword arguments:
    turbine -- turbine object
    point -- point to coat
    """

    with turbine.robot:
        iksol, tolerance = ik_angle_tolerance_normal_plane(turbine, point, iter_surface)
        feasible_iksol = []
        feasible_tolerance = []
        for i in range(0,len(iksol)):
            if(trajectory_constraints(turbine, iksol[i], point)):
                feasible_iksol.append(iksol[i])
                feasible_tolerance.append(tolerance[i])
        return feasible_iksol, feasible_tolerance

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
    if not turbine.env.GetCollisionChecker().SetCollisionOptions(CollisionOptions.Distance):
        collisionChecker = RaveCreateCollisionChecker(turbine.env,'pqp')
        collisionChecker.SetCollisionOptions(CollisionOptions.Distance)
        turbine.env.SetCollisionChecker(collisionChecker)
        
    robot = turbine.robot
    q0 = robot.GetDOFValues()
    manip = robot.GetActiveManipulator()
    lower_limits, upper_limits = robot.GetActiveDOFLimits()
    lower_limits = lower_limits+0.01
    upper_limits = upper_limits-0.01
    report = CollisionReport()
    
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
        
        def self_collision_cons(q):
            robot.SetDOFValues(q)
            report = CollisionReport()
            turbine.robot.CheckSelfCollision(report)
            return report.minDistance-0.007
        
        def env_collision_cons(q):
            robot.SetDOFValues(q)
            report = CollisionReport()
            turbine.env.CheckCollision(turbine.robot,report)
            return report.minDistance-0.007

        
            ## New opt
            #pn1 = mathtools.compute_perpendicular_vector(point[3:6])
            #pn2 = cross(pn1,point[3:6])
            #pn3 = point[3:6]
            #pos_on_pn3 = dot(pos,point[3:6])*pos
            #if linalg.norm(pos_on_pn3 - point[0:3])<=0.03:
            #    pos_on_pn3 = point[0:3]
            #pos = dot(pos,pn1)*pos + dot(pos,pn2)*pos + dot(pos,pn3)*pos
            #dif = pos-point[0:3]
            #return dot(dif,dif)/tol            
        
        cons = ({'type':'eq', 'fun': position_cons},
                {'type':'ineq', 'fun': self_collision_cons},
                {'type':'ineq', 'fun': env_collision_cons})
        
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

    for angle in arange(0, tolerance, angle_discretization):
        normal_tol = dot(point[3:6],transpose(mathtools.Raxis(tan, angle)))
        sols = ikfast(robot, concatenate((point[0:3], normal_tol)))
        for sol in sols:
            iksol.append(sol)
            angle_tolerance.append(angle)
            
        normal_tol = dot(point[3:6],transpose(mathtools.Raxis(tan, -angle)))
        sols = ikfast(robot, concatenate((point[0:3], normal_tol)))
        for sol in sols:
            iksol.append(sol)
            angle_tolerance.append(angle)
            
    return iksol, angle_tolerance 


def ikfast(robot, point):
    """
    Call openrave IKFast. It computes the inverse kinematic for the point.
    It returns all solutions.

    keyword arguments:
    robot -- the robot. 
    point -- point to coat is a 6D array, which (x,y,z) cartesian position
    and (nx,ny,nz) is the normal vector of the point, w.r.t. the world frame.
    """
    point = array(point)
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

def compute_manipulability_det(joint_configuration, robot):
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
    return sqrt(linalg.det(dot(transpose(J),J)))#, sqrt(linalg.det(dot(Jpos,transpose(Jpos)))), sqrt(linalg.det(dot(Jori,transpose(Jori))))

def best_joint_solution_regarding_manipulability(joint_solutions, tolerance, robot):
    """
    Given a list of joint solutions for a specific point, the function computes
    manipulability and returns the best solution regarding manipulability criteria.
    """
                                 
    biggest_manipulability = 0
    temp_q = []
    tolerance = array(tolerance)
    joint_solutions = array(joint_solutions)
    joint_solutions = joint_solutions[abs(tolerance)==min(abs(tolerance))]
    
    for q in joint_solutions:
        try:
            manipulability = compute_manipulability_det(q, robot)
            if manipulability > biggest_manipulability:
                biggest_manipulability = manipulability
                temp_q = q
        except IndexError:
            raise IndexError('There is no solution for the given trajectory.')
    return temp_q

def trajectory_constraints(turbine, iksol, point):
    """
    Check robot self and environment collision, optimization result,
    and angle tolerance.
    """

    robot = turbine.robot

    # Verifying environment and self collision
    with robot:
        robot.SetDOFValues(iksol)
        
        # Verifying angle tolerance
        if not orientation_cons(turbine, point):
            return False
    return True

def set_ikmodel_transform6D(robot):
    ikmodel = databases.inversekinematics.InverseKinematicsModel(
        robot=robot, iktype=IkParameterization.Type.Transform6D)
    if not ikmodel.load():
        ikmodel.autogenerate()
    return  

def random_joint(robot, joint):
    rjoints = joint+random.rand(len(robot.GetDOFValues()))
    limitsdow = robot.GetActiveDOFLimits()[0]
    limitsup = robot.GetActiveDOFLimits()[1]
    for j in range(0,len(rjoints)):
        rjoints[j] = min([rjoints[j],limitsup[j]])
        rjoints[j] = max([rjoints[j],limitsdow[j]])
    return rjoints

def generate_random_joint_solutions(turbine, point, tries):
    robot = turbine.robot
    joints = []
    for i in range(0,tries):
        random_joint(robot, zeros(len(robot.GetDOFValues())))
        res = orientation_error_optimization(turbine, point)
        if res.success:
            if trajectory_constraints(turbine, res.x, point):
                break
    else:
        return joints

    jointx = res.x
    joints.append(jointx)
    with robot:
        for i in range(0,tries):
            rjoints = random_joint(robot, jointx)
            robot.SetDOFValues(rjoints)
            res = orientation_error_optimization(turbine, point)
            if res.success:
                if trajectory_constraints(turbine, res.x, point):
                    if list(res.x) not in [list(x) for x in joints]:
                        joints.append(res.x)
    return joints
       
def joint_distance(joint1, joint2):
    dif = abs(joint1-joint2)
    return sum(dif)

def joint_planning(turbine, ordered_waypoints, tries = 12):
    joints = []
    robot = turbine.robot
    for i in range(0,len(ordered_waypoints)):
        joints.append(ikfast(robot, ordered_waypoints[i]))

    for i in range(0,len(joints)):
        if len(joints[i])==0:
            if i==0:
                joints[i] = generate_random_joint_solutions(turbine, ordered_waypoints[i], tries)
            else:
                joints_sol = []
                for joint in joints[i-1]:
                    robot.SetDOFValues(joint)
                    res = orientation_error_optimization(turbine, ordered_waypoints[i])
                    if res.success:
                        if trajectory_constraints(turbine, res.x, ordered_waypoints[i]):
                            joints_sol.append(res.x)
                joints[i] = joints_sol

    for i in range(0,len(joints)):
        if len(joints[i])==0:
            print ordered_waypoints[i]
            raise IndexError('joints with zero length')
    
    return joints

def compute_foward_cost(joints0,joints1):
    cost = zeros((len(joints0),len(joints1)))
    for i in range(0,len(joints0)):
        for j in range(0,len(joints1)):
            cost[i][j] = joint_distance(joints0[i], joints1[j])
    return cost

def make_dijkstra(joints, verbose = False):
    virtual_start = (-1,-1)
    virtual_end = (-2,-2)
    adj = dict()
    cost = dict()
    for jointsi in range(0,len(joints)-1):
         mcost = compute_foward_cost(joints[jointsi],joints[jointsi+1])
         for u in range(0,len(mcost)):
             for v in range(0,len(mcost[u])):
                 l = adj.get((jointsi,u),[])
                 l.append((jointsi+1,v))
                 adj[(jointsi,u)] = l
                 cost[((jointsi,u),(jointsi+1,v))] = mcost[u][v]

    for joints0i in range(0,len(joints[0])):
        l = adj.get(virtual_start,[])
        l.append((0,joints0i))
        adj[virtual_start] = l
        cost[(virtual_start,(0,joints0i))] = 0

    for jointsi in range(0,len(joints[-1])):
        l = adj.get(virtual_end,[])
        l.append((len(joints)-1,jointsi))
        adj[virtual_end] = l
        cost[(virtual_end,(len(joints)-1,jointsi))] = 0

        l = adj.get((len(joints)-1,jointsi),[])
        l.append(virtual_end)
        adj[(len(joints)-1,jointsi)] = l
        cost[((len(joints)-1,jointsi),virtual_end)] = 0

    cost = dijkstra2.make_undirected(cost)
    predecessors, min_cost = dijkstra2.dijkstra(adj, cost, virtual_start, virtual_end)
    c = virtual_end
    path = [c]
    
    while predecessors.get(c):
        path.insert(0, predecessors[c])
        c = predecessors[c]

    joint_path = []
    for i in range(1,len(path)-1):
        joint_index = path[i][0]
        joint_configuration = path[i][1]
        joint_path.append(joints[joint_index][joint_configuration])

    if verbose:
        return joint_path, path, min_cost, adj, cost
    else:
        return joint_path






