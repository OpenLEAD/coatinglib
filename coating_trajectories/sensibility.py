#!/usr/bin/env python
import db
from os.path import join, realpath
from os import makedirs, environ
import blade_modeling
from turbine import Turbine
from turbine_config import TurbineConfig, ConfigFileError
import planning
from visualizer import Visualizer
from numpy import all, nonzero, split, array
import cPickle
import errno
import rail_place
from itertools import izip

def get_joint(joint_base, point_num):
        return joint_base[point_num]

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
            w_list, alpha_list, times = planning.compute_general_velocities(turb, joints[i][j], trajectories[i][j])
            torque = []
            sensible = []
            for ray, joints, w, alpha in izip(trajectories[i][j], joints[i][j], w_list, alpha_list):
                torque += [planning.torque_computation(turb, joints, w, alpha)]
                sensible += [planning.sensibility(turb, ray, w, alpha)]
                
            velocity_times += [(w_list, alpha_list, times)]
            torque_seg += [torque]
            sensible_seg += [sensible]
            
        velocities_times += [velocity_times]
        torque_traj += [torque_seg]
        sensible_traj += [sensible_seg]
    return velocities_times, torque_traj, sensible_traj

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
            w_list, alpha_list, times = planning.compute_general_velocities(turb, joints[i][j], trajectories[i][j])
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

def plot_segs_comp(base_num):
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
