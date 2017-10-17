#!/usr/bin/env python
import db
from os.path import join
from os import makedirs
from numpy import all, nonzero, split, mean, cumsum
from numpy import sqrt, dot, concatenate, array, cross
from numpy import abs, minimum, arccos
from scipy.optimize import linprog
from numpy.linalg import norm
from mathtools import MLS, roundrobin, general_finite_difference, hat
import errno
import robot_utils
import rail_place
from itertools import izip

def get_joint(joint_base, point_num):
        return joint_base[point_num]


def MLS_general_velocities(turbine, joints_trajectory, points, n=4, scale=1.5):
    """
    points = rays[:,0:3]
    0)Compute times based on h and distances calculated from rays

    1)Use MLS analytical derivative to estimate velocities and accelerations.

    """
    h = turbine.config.coating.coating_speed
    dtimes = norm(array(points[1:]) - array(points[:-1]), axis=1) / h
    scale *= mean(dtimes)
    times = cumsum(array([0.] + list(dtimes)))

    w_list = []
    alpha_list = []

    for joints in array(joints_trajectory).T:
        _, w, alpha = MLS(joints, times, times, n, scale)
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
    dpoints = array(points[1:]) - array(points[:-1])
    times = cumsum(array([0.] + list(norm(dpoints, axis=1))) / h)

    N = len(times)

    w_list = []
    alpha_list = []

    for i, time in enumerate(times):
        past = range(i)[::-1]
        future = range(i, N)
        order = list(roundrobin(future, past))
        order = order[:min(len(order), 7)]
        _, w, alpha = general_finite_difference(time, array(joints_trajectory)[order], times[order])
        w_list += [w]
        alpha_list += [alpha]
    return w_list, alpha_list, times

def check_velocities_DB(DB, base):
    segments = db.load_pickle(join(DB.db_main_path,'seg',str(base)+'.pkl'))[base]
    joints_DB = db.load_pickle(join(DB.db_main_path,'joints',str(base)+'.pkl'))[base]
    blade = DB.load_blade()
    ntp = DB.get_sorted_points()
    trajectories = []
    joints = []
    for segment in segments:
        trajectory = []
        joint = []
        for seg in segment:
            trajectory.append(map( lambda x: blade.compute_ray_from_point(ntp[x]), seg ))
            joint.append(map( lambda x: get_joint(joints_DB,x), seg ))
        trajectories.append(trajectory)
        joints.append(joint)    
    velocities = compute_velocities(trajectories, joints, DB.turb) 
    return check_angular_velocities_segs(velocities, segments, joints, base, DB.turb)
  
def compute_sensibility(trajectories, joints, turb):
    """
    trajectories -- Trajectories of the segments Ni x Mj x Ok x 6
    Ni (parallels) x Mj (segments inside each parallel) x Ok (rays inside each segment) x 6 (ray size):
    [  [[],[],[], ...        ], [[],[], ...            ],...]
      |segments of parallel 0|, |segments of parallel 1| 
    joints -- joints Ni x Mj x Ok x 6 according to segments
    """
    velocities_times = []
    torque_traj = []
    sensible_traj = []

    for i in range(0,len(trajectories)):
        velocity_times = []
        torque_seg = []
        sensible_seg = []
        for j in range(0,len(trajectories[i])):
            w_list, alpha_list, times = compute_general_velocities(turb, joints[i][j], trajectories[i][j])
            torque = []
            sensible = []
            for ray, joints, w, alpha in izip(trajectories[i][j], joints[i][j], w_list, alpha_list):
                torque += [torque_computation(turb, joints, w, alpha)]
                sensible += [sensibility(turb, ray, w, alpha)]
                
            velocity_times += [(w_list, alpha_list, times)]
            torque_seg += [torque]
            sensible_seg += [sensible]
            
        velocities_times += [velocity_times]
        torque_traj += [torque_seg]
        sensible_traj += [sensible_seg]
    return velocities_times, torque_traj, sensible_traj


def torque_computation(turbine, joints, w, alpha, verify=False):
    """
    Apply joints position 'joints', joints velocities 'w' and joints acceleration 'alpha'
    Return torques - (6,) array
    verify - if True function returns False when a torque is bigger than maximum
    """
    with turbine.robot:
        gun = turbine.robot.GetLink('Gun')
        flame_force = {
            gun.GetIndex(): list(-gun.GetTransform()[0:3, 2] * turbine.config.coating.flame_thrust) + [0, 0, 0]}
        turbine.robot.SetDOFValues(joints, turbine.robot.GetActiveDOFIndices(),
                                   turbine.robot.CheckLimitsAction.CheckLimitsSilent)
        turbine.robot.SetDOFVelocities(w, turbine.robot.CheckLimitsAction.CheckLimitsSilent,
                                       turbine.robot.GetActiveDOFIndices())
        torques = turbine.robot.ComputeInverseDynamics(alpha, flame_force)

    if any(torques > turbine.robot.GetDOFMaxTorque) and verify:
        return False
    else:
        return torques


def sensibility(turbine, ray, w, alpha):
    normal_vec = ray[3:6]
    tangent_vec = cross(ray[0:3], normal_vec)
    tangent_vec = tangent_vec / sqrt(dot(tangent_vec, tangent_vec))
    perp_vec = cross(normal_vec, tangent_vec)

    Jpos = turbine.manipulator.CalculateJacobian()
    Hpos = turbine.robot.ComputeHessianTranslation(turbine.manipulator.GetArmDOF(),
                                                   turbine.manipulator.GetEndEffectorTransform()[0:3, 3])
    Hpos = dot(Hpos, w)
    theta_limits = zip(-turbine.robot.GetDOFResolutions(), turbine.robot.GetDOFResolutions())
    w_limits = zip(-abs(w) * 0.001, abs(w) * 0.001)  # HARDCODED 1% error
    limits = tuple(theta_limits + w_limits)

    Hpos_tan = dot(Hpos, tangent_vec)
    Jpos_tan = dot(tangent_vec, Jpos)
    errorgain_tan = concatenate((Jpos_tan, Hpos_tan))
    velocity_tan_error = (
    linprog(errorgain_tan, bounds=limits).get('fun'), -linprog(-errorgain_tan, bounds=limits).get('fun'))

    Jpos_normal = dot(normal_vec, Jpos)
    position_normal_error = (
    linprog(Jpos_normal, bounds=theta_limits).get('fun'), -linprog(-Jpos_normal, bounds=theta_limits).get('fun'))

    Jpos_perp = dot(perp_vec, Jpos)
    position_perp_error = (
    linprog(Jpos_perp, bounds=theta_limits).get('fun'), -linprog(-Jpos_perp, bounds=theta_limits).get('fun'))

    Jw = turbine.manipulator.CalculateAngularVelocityJacobian()
    x_dir = turbine.manipulator.GetTransform()[0:3, 0]
    nhat = hat(normal_vec)
    xn = dot(x_dir, nhat)
    Jcos = -dot(xn, Jw)
    cosn = -dot(x_dir, normal_vec)
    dcosn = (linprog(Jcos, bounds=theta_limits).get('fun'), -linprog(-Jcos, bounds=theta_limits).get('fun'))
    angle_error = tuple(arccos(minimum(cosn + dcosn, 1.)) - arccos(cosn))

    return velocity_tan_error, position_normal_error, position_perp_error, angle_error

def compute_velocities(trajectories, joints, turb):
    """
    trajectories -- Trajectories of the segments Ni x Mj x Ok x 6
    Ni (parallels) x Mj (segments inside each parallel) x Ok (rays inside each segment) x 6 (ray size):
    [  [[],[],[], ...        ], [[],[], ...            ],...]
      |segments of parallel 0|, |segments of parallel 1| 
    joints -- joints Ni x Mj x Ok x 6 according to segments
    """
    velocities_times = []
    torque_traj = []
    sensible_traj = []

    for i in range(0,len(trajectories)):
        velocity_times = []
        for j in range(0,len(trajectories[i])):
            w_list, alpha_list, times = compute_general_velocities(turb, joints[i][j], trajectories[i][j])
            velocity_times += [(w_list, alpha_list, times)]          
        velocities_times += [velocity_times]
    return velocities_times


def check_angular_velocities(velocities, turb):
    maxVel = turb.robot.GetDOFVelocityLimits()
    return all((abs(velocities)-maxVel)<0,1)

def check_angular_velocities_segs(velocity_segs, segments, joints, base_num, turb):

    new_seg_base = dict()
    new_seg_base[base_num] = list()
    
    for i,velocity_seg in enumerate(velocity_segs):
        if len(velocity_seg)>0:
            new_seg = []
            for j,vel_alpha in enumerate(velocity_seg):
                vel = vel_alpha[0]
                check = check_angular_velocities(array(vel), turb)
                indices = nonzero(check[1:] != check[:-1])[0] + 1
                b = split(array(segments[i][j]), indices)
                b = b[0::2] if check[0] else b[1::2]
                for bi in b:
                    if len(bi)>=3:
                        new_seg.append(list(bi))
            if len(new_seg)>0:
                new_seg_base[base_num].append(new_seg)
    return new_seg_base

def compute_new_segs_DB(DB):
    new_path_seg = join(DB.db_main_path,'new_seg')

    try:
        makedirs(new_path_seg)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
    
    while True:
        visited_bases = db.load_pickle(join(DB.db_main_path,'visited_bases.pkl'))
        base_num = None
        for key, value in visited_bases.iteritems():
            if value == False:
                visited_bases[key] = True
                base_num = key
                break
        if base_num is None:
            break
        db.save_pickle(visited_bases,join(DB.db_main_path,'visited_bases.pkl'))
        del visited_bases

        try:
            seg_base = db.load_pickle(join(DB.db_main_path,'seg',str(base_num)+'.pkl'))
        except IOError:
            continue

        new_seg_base = check_velocities_DB(DB, base_num)

        print 'saving base_num: ', base_num

        db.save_pickle(new_seg_base, join(new_path_seg,str(base_num)+'.pkl'))
    return


def best_joint_solution_regarding_manipulability(joint_solutions, tolerance, robot):
    """
    Given a list of joint solutions for a specific point, the function computes
    manipulability and returns the best solution regarding manipulability criteria.
    """

    biggest_manipulability = 0
    temp_q = []
    tolerance = array(tolerance)
    joint_solutions = array(joint_solutions)
    joint_solutions = joint_solutions[abs(tolerance) == min(abs(tolerance))]

    for q in joint_solutions:
        try:
            manipulability = robot_utils.compute_manipulability_det(q, robot)
            if manipulability > biggest_manipulability:
                biggest_manipulability = manipulability
                temp_q = q
        except IndexError:
            raise IndexError('There is no solution for the given trajectory.')
    return temp_q


def plot_segs_comp(DB,vis,turb,base_num,segs_path,new_segs_path):
    seg_base = DB.load_db_pickle(join(segs_path,str(base_num)+'.pkl'))
    new_seg_base = DB.load_db_pickle(join(new_segs_path,str(base_num)+'.pkl'))
    trajectories = seg_base[base_num]
    new_trajectories = new_seg_base[base_num]
    ntp = DB.get_sorted_points()

    for i in range(0,len(trajectories)):
        if len(trajectories[i])>0:
            for j in range(0,len(trajectories[i])):
                for num_point in trajectories[i][j]:
                    p = vis.plot(ntp[num_point],'p',(1,0,0))

    for i in range(0,len(new_trajectories)):
        if len(new_trajectories[i])>0:
            for j in range(0,len(new_trajectories[i])):
                for num_point in new_trajectories[i][j]:
                    p = vis.plot(ntp[num_point],'p',(0,0,1))

    base = DB.get_sorted_bases()[base_num]
    rp = rail_place.RailPlace(base)
    turb.place_rail(rp)
    turb.place_robot(rp)
    return 
