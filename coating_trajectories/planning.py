from numpy import sqrt, dot, concatenate, array, transpose, linalg, cross, zeros, eye, max
from numpy import abs, cumsum, minimum, arccos, random, linspace, mean, inf
from numpy.linalg import norm
from openravepy import IkFilterOptions, interfaces, databases, IkParameterization
from openravepy import CollisionOptions, RaveCreateCollisionChecker, CollisionReport
from math import pi, cos, sin
from scipy.optimize import minimize, linprog
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

def MLS_general_velocities(turbine, joints_trajectory, points, n = 4, scale = 1.5):
    """
    points = rays[:,0:3]
    0)Compute times based on h and distances calculated from rays

    1)Use MLS analytical derivative to estimate velocities and accelerations.
    
    """
    h = turbine.config.coating.coating_speed
    dtimes = norm(array(points[1:])-array(points[:-1]),axis=1)/h
    scale *= mean(dtimes)
    times = cumsum(array([0.]+list(dtimes)))

    w_list = []
    alpha_list = []

    
    for joints in array(joints_trajectory).T:
        _,w,alpha = mathtools.MLS(joints, times,times,n,scale)
        w_list += [w]
        alpha_list += [alpha]
    
    
    return list(array(w_list).T), list(array(alpha_list).T), times
    
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
    
def compute_robot_joints(turbine, trajectory, deep=False):
    """
    Iterates points of the trajectories, computing optimal robot's joints
    (minimizing orientation error). For the first point, it uses ikfast.

    Args:
        turbine: (@ref Turbine) turbine object.
        trajectory: trajectory to coat
        trajectory_index: where to begin. Index of a feasible point in the trajectory.
    """

    joint_solutions = []

    # Find solutions for points
    for index in range(0, len(trajectory)):
        iksol = ik_angle_tolerance(turbine, trajectory[index], deep=deep)
        if len(iksol) == 0:
            return joint_solutions
        else:
            joint_solutions.append(iksol)
    return joint_solutions


def compute_first_feasible_point(turbine, trajectory):
    """ Method to compute the first feasible point in the trajectory: where to start.

    Args:
        turbine: (@ref Turbine) turbine object.
        trajectory: trajectory to coat
    """

    robot = turbine.robot
    with robot:
        for i in range(0,len(trajectory)):
            sols = ik_angle_tolerance(turbine, trajectory[i], deep=False)
            if len(sols)>0:
               return i, sols
        raise ValueError('No solution for given trajectory')

def orientation_cons(turbine, point):
    """
    Return True if the orientation error is inside limits.
    """

    Rx = turbine.robot.GetActiveManipulator().GetTransform()[0:3,0]
    Rx = Rx/linalg.norm(Rx)
    return (dot(point[3:6], Rx)) + cos(turbine.config.coating.angle_tolerance) <= 0


def ik_angle_tolerance(turbine, point, angle_tolerance_init=0, angle_tolerance_end=None, number_of_phi = 24, number_of_theta=7, deep=False):
    """ Solve the inverse kinematics given point (IKFast) with maximum tolerance angle.

    Args:
        turbine: (@ref Turbine) turbine object.
        point: point to be coated is a 6D array, which (x,y,z) cartesian position
        and (nx,ny,nz) is the normal vector of the point, w.r.t. the world frame.
    """

    robot = turbine.robot
    iksol = []
    
    # Compute solution with maximum distance and angle tolerances
    if angle_tolerance_end is None: 
        angle_tolerance_end = turbine.config.coating.angle_tolerance

    for theta in linspace(angle_tolerance_init, angle_tolerance_end, number_of_theta):
        normal_tol = dot(point[3:6],transpose(mathtools.Raxis(
            mathtools.compute_perpendicular_vector(point[3:6]), theta)))
        normal_tol = normal_tol/linalg.norm(normal_tol)
        for phi in linspace(0,2*pi,number_of_phi*sin(theta)+1,endpoint=False):
            iksoli = ikfast(robot, concatenate((point[0:3],dot(
                normal_tol, transpose(mathtools.Raxis(point[3:6],phi))))))
            iksol.extend(iksoli)
        if not deep:
            if len(iksol)>0:
                return iksol
    return iksol

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
        solutions = robot.GetActiveManipulator().FindIKSolutions(T, True)

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

def single_vel_distance( vel1, vel2, dt, vel_limits, acc_limits):
    dif = abs(array(vel1)-array(vel2))/dt
    return sum((dif/acc_limits)**2) + sum((abs(vel2)/vel_limits)**2)
       
def joint_distance_mh12(joint1, joint2, dt, vel_limits):
    dif = abs(array(joint1) - array(joint2))
    if dt == 0:
        dif = max(dif,1)
        blown = (dif > 1e-5)
        dif[blown] = inf
        dif[~blown] = 0
        return dif
    dif /= dt
    percent_dif = dif/vel_limits
    return max(percent_dif,1)

def compute_foward_cost(joints0, joints1, limits):
    cost = zeros((len(joints0),len(joints1)))
    for i in range(0,len(joints0)):
        cost[i] = joint_distance_mh12(joints0[i], joints1, limits)
    return cost


def make_dijkstra(joints, dtimes, vel_limits, acc_limits, verbose = False):
    virtual_start = (-1,-1)
    virtual_end = (-2,-1)
    adj = dijkstra2.dijkstra_adj(joints,dtimes)

    for jointsi in range(len(joints)-1):
        for u in range(len(joints[jointsi])):
            for v in range(len(joints[jointsi+1])):
                adj.add_link((jointsi,u),(jointsi+1,v))

    for joints0i in range(len(joints[0])):
        adj.add_link(virtual_start,(0,joints0i))

    for jointsi in range(len(joints[-1])):
        adj.add_link((len(joints)-1,jointsi),virtual_end)

    vs = dijkstra2.virtual_node((-1,-1),tuple(zeros(len(joints[0][0]))))
    ve = dijkstra2.virtual_node((-2,-1),tuple(zeros(len(joints[0][0]))))
    dtimes[0] = 1

    cost = dijkstra2.dijkstra_acc_cost(single_vel_distance,vs,ve,dtimes,vel_limits,acc_limits)

    predecessors, min_cost = dijkstra2.dijkstra(adj, cost, vs, ve)

    c = next(y for y in predecessors.keys() if y == ve)

    path = [c]
    while predecessors.get(c):
        path.insert(0, predecessors[c])
        c = predecessors[c]

    joint_path = []
    for i in range(1,len(path)-1):
        joint_index = path[i][0][0]
        joint_configuration = path[i][0][1]
        joint_path.append(joints[joint_index][joint_configuration])

    if verbose:
        return joint_path, path, min_cost, adj, cost
    else:
        return joint_path


def make_dijkstra_vel(joints, dtimes, vel_limits, acc_limits, verbose = False):
    virtual_start = (-1,-1)
    virtual_end = (-2,-2)
    adj = dict()

    for jointsi in range(0,len(joints)-1):
         for u in range(0,len(joints[jointsi])):
             for v in range(0,len(joints[jointsi+1])):
                 l = adj.get((jointsi,u),[])
                 l.append((jointsi+1,v))
                 adj[(jointsi,u)] = l

    for joints0i in range(0,len(joints[0])):
        l = adj.get(virtual_start,[])
        l.append((0,joints0i))
        adj[virtual_start] = l

    for jointsi in range(0,len(joints[-1])):
        l = adj.get((len(joints)-1,jointsi),[])
        l.append(virtual_end)
        adj[(len(joints)-1,jointsi)] = l

    cost = dijkstra2.dijkstra_vel_cost(joint_distance_mh12,joints,virtual_start,virtual_end,dtimes,vel_limits)

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