from numpy import sqrt, dot, concatenate, arange, array, zeros, transpose, linalg
from openravepy import IkFilterOptions
from math import pi, cos, sin, atan2
from scipy.optimize import minimize
from copy import copy
from collections import deque
import mathtools
import time

""" Main package for robot joints' positions and velocities planning,
robot base calculation, torque and manipulability analysis.

Keyword arguments:
places - robot places to coat a specific set of parallels
turbine - full environment and parameters
"""

_filter_options = {'mh12': mh12_filter}

def _std_robot_filter(turbine, trajectories):
    raise ValueError("No trajectory filter for "+turbine.robot.GetName()+" robot. Create new function.")

def mh12_filter(turbine, trajectories):
    _working_radius_squared = 1.285**2
    def distance_robot_squared(point,turbine):
        delta = (tubine.robot.GetJoints[1].GetAnchor()-point)
        return dot(delta,delta)
    
    filtered_trajectories = []
    for trajectory in trajectories:
        
        filtered_trajectory = []
        for point in trajectory:
            
            if distance_robot_squared(point[0:3],turbine) < _working_radius_squared:
                filtered_trajectory += [point]
        if not filtered_trajectory:
            filtered_trajectories += [filtered_trajectory]

    return filtered_trajectories
    

def filter_trajectories(turbine, trajectories):
    name = turbine.robot.GetName()
    _filter_options.get(name,_std_robot_filter)(turbine, trajectories)

def workspace_limit(turbine,trajectory, joints_trajectory, trajectory_index):

    # 'j' for joints, so it doesnt clumpsy the equations
    j = deque(joints_trajectory)
    j.rotate(-trajectory_index)

    h = turbine.config.model.trajectory_step
    
    # Joints velocity - Central Difference (h**6 order error)
    w = ( (j[3]-j[-3]) + 9*(j[-2]-j[2]) + 45*(j[1]-j[1]) )/(60.0*h)

    # Joints acceleration - Central Difference (h**6 order error)
    alpha = ( 2*(j[-3]+j[3]) - 27*(j[-2]+j[2]) + 270*(j[-1]+j[1]) - 490*j[0] )/(180.0*h**2)

    # INVERSE DYNAMICS ComputeInverseDynamics
    with turbine.robot:
        turbine.robot.SetDOFValues(j[0], _, turb.robot.CheckLimitsAction.CheckLimitsThrow)
        turbine.robot.SetDOFVelocities(w, _, turb.robot.CheckLimitsAction.CheckLimitsThrow)
        torques = turb.robot.ComputeInverseDynamics(alpha)
        if all(torques < turb.robot.GetDOFMaxTorque):
            return False#THROW EXCEPTION
        
    

def compute_robot_joints(turbine, trajectory):
    """
    Iterates points of the trajectories, computing optimal robot's joints
    (minimizing orientation error).

    Keyword arguments:
    place -- robot position
    turbine -- tubine object
    """

    robot = turbine.robot
    q0, _ = inverse_kinematics(turbine, trajectory[0])
    biggest_manipulability = 0
    temp_q = []
    for q in q0:
        try: 
            manipulability, _, _ = compute_manipulability_det(robot, q)
            if manipulability>biggest_manipulability:
                biggest_manipulability = manipulability
                temp_q = q
        except IndexError:
            raise IndexError('There is no solution for the given trajectory.')
    
    robot.SetDOFValues(temp_q)
    joint_solutions = []
    
    for index in range(0,len(trajectory)):
        res = orientation_error_optimization(turbine, trajectory[index])
        if not res.success:
            print 'Trajectory terminates in point: '+str(index)+'/'+str(len(trajectory))
            return joint_solutions, False 
        else:
            joint_solutions.append(res.x)
            robot.SetDOFValues(res.x)
            time.sleep(0.01)
    return joint_solutions, True

def orientation_error_optimization(turbine, point):
    """ Minimizes the orientation error, given an initial configuration of the robot
    and the desired point to coat (with its normal vector). The constraints of the minimi-
    zation are: to reach de point position (equality), to be in coating angle tolerance
    (inequality).

    Keyword arguments:
    robot -- the robot.
    point -- 6x1 vector is the goal (point to be coated). (x,y,z) cartesian position
    and (nx,ny,nz) is the normal vector of the point, w.r.t. the world frame.
    """
    #TODO check collision
    robot = turbine.robot
    cos_tolerance = cos(turbine.config.coating.angle_tolerance)
    q0 = robot.GetDOFValues()
    manip = robot.GetActiveManipulator()
    lower_limits, upper_limits = robot.GetActiveDOFLimits()
    lower_limits = lower_limits+0.05
    upper_limits = upper_limits-0.05
    
    with robot:
        def func(q):
            robot.SetDOFValues(q, manip.GetArmIndices())
            T = manip.GetTransform()
            Rx = dot(manip.GetDirection(), T[0:3,0:3])
            Rx = Rx/sqrt(dot(Rx,Rx))
            return -dot(point[3:6], Rx)
        def position_cons(q):
            robot.SetDOFValues(q, manip.GetArmIndices())
            T = manip.GetTransform()
            pos = T[0:3,3]
            dif = pos-point[0:3]
            return dot(dif,dif)
        def orientation_cons(q):
            robot.SetDOFValues(q, manip.GetArmIndices())
            T = manip.GetTransform()
            Rx = dot(manip.GetDirection(), T[0:3,0:3])
            Rx = Rx/sqrt(dot(Rx,Rx))
            return dot(point[3:6], Rx) - cos_tolerance
        cons = ({'type':'eq', 'fun': position_cons},
                {'type':'ineq', 'fun': orientation_cons})
        bnds = tuple([(lower_limits[i],upper_limits[i]) for i in range(0,len(lower_limits))])
        res = minimize(func, q0, constraints=cons, method='SLSQP',
                       bounds = bnds, options={'disp': False})
    return res    

def Qvector_backtrack(y, Q):
    """ Call when optimization fails. Previous solution should be recomputed.

    Keyword arguments:
    y -- computed coated points.
    Q -- computed robot's joints.
    """
    for rev_y in reversed(y):
        res = coating.optmizeQ(turbine, rev_y, Q[-1])
        Q.append(res.x)
        if not res.success:
            print 'Qvector_backtrack error'
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
    if len(iksol)>0:return iksol, []

    # Compute solution with distance and angle tolerances
    angles = []
    if turbine.config.coating.angle_tolerance>0:
        angles = arange(0, turbine.config.coating.angle_tolerance, 0.001)
    number_of_angles = 10

    if (coating_tolerance_max!=0) or (coating_tolerance_min!=0):
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
        return robot.GetActiveManipulator().FindIKSolutions(T, IkFilterOptions.CheckEnvCollisions)

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
