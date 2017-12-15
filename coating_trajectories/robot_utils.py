from numpy import sqrt, dot, concatenate, array, transpose, linalg, cross, zeros, eye, max
from numpy import abs, minimum, arccos, random, linspace, inf
from openravepy import IkFilterOptions, interfaces, databases, IkParameterization
from openravepy import CollisionOptions, RaveCreateCollisionChecker, CollisionReport
from math import pi, sin
from mathtools import central_difference
import mathtools

"""
Main package for robot joints' positions and velocities planning,
robot base calculation, torque and manipulability analysis.
"""


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
        Rx = Tee[0:3, 0]
        Rx = Rx / linalg.norm(Rx)
        Rab = mathtools.Rab(Rx, -point[3:6])

        T = eye(4)
        T[0:3, 0:3] = dot(Rab, Tee[0:3, 0:3])
        T[0:3, 3] = point[0:3]
        solutions = robot.GetActiveManipulator().FindIKSolutions(T, True)

        if len(solutions) > 0:
            if len(solutions.shape) == 1:
                solutions = solutions.reshape((1, solutions.shape[0]))

        return solutions

def compute_angular_velocities(turbine, joints_trajectory, trajectory_index):

    if (trajectory_index>2) and ((len(joints_trajectory)-trajectory_index)>3):
        return central_difference(turbine, joints_trajectory, trajectory_index)
    else: return None


def ik_angle_tolerance(turbine, point, angle_tolerance_init=0, angle_tolerance_end=None, number_of_phi=24,
                       number_of_theta=7, deep=False):
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
        normal_tol = dot(point[3:6], transpose(mathtools.Raxis(
            mathtools.compute_perpendicular_vector(point[3:6]), theta)))
        normal_tol = normal_tol / linalg.norm(normal_tol)
        for phi in linspace(0, 2 * pi, number_of_phi * sin(theta) + 1, endpoint=False):
            iksoli = ikfast(robot, concatenate((point[0:3], dot(
                normal_tol, transpose(mathtools.Raxis(point[3:6], phi))))))
            iksol.extend(iksoli)
        if not deep:
            if len(iksol) > 0:
                return iksol
    return iksol


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
        J = concatenate((Jpos, Jori))
    return sqrt(linalg.det(dot(transpose(J),
                               J)))  # , sqrt(linalg.det(dot(Jpos,transpose(Jpos)))), sqrt(linalg.det(dot(Jori,transpose(Jori))))


def compute_robot_joints(turbine, trajectory, deep=False):
    """ Iterates points of the trajectories. It uses ikfast.

    Args:
        turbine: (@ref Turbine) turbine object.
        trajectory: trajectory to coat
        trajectory_index: where to begin. Index of a feasible point in the trajectory.
    """

    joint_solutions = []

    # Find solutions for points
    for index in range(0, len(trajectory)):
        iksol = ik_angle_tolerance(turbine, trajectory[index], deep=deep)
        if iksol:
            joint_solutions.append(iksol)
        else:
            return None
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