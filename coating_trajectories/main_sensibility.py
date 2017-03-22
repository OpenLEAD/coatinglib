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

def load_blade(folder):
    """
    Function to load blade model given folder.
    """
    xml_trajectories_path = join(folder,"trajectory/trajectory.xml")
    blade = blade_modeling.BladeModeling(turb, turb.blades[0])
    blade.load_trajectory(xml_trajectories_path)
    return blade

def get_joint(base_num,joint_base, point_num):
        return joint_base[base_num][point_num]
    
def compute_velocities(base_num,turbine):
    blade = load_blade(blade_folder)
    
    seg_base = DB.load_db_pickle(join(segs_path,str(base_num)+'.pkl'))
    joint_base = DB.load_db_pickle(join(joints_path,str(base_num)+'.pkl'))
    trajectories = seg_base[base_num]
    velocities_times = []
    torque_traj = []
    sensible_traj = []

    for i in range(0,len(trajectories)):
        velocity_times = []
        torque_seg = []
        sensible_seg = []
        if len(trajectories[i])>0:
            list_rays = DB.compute_rays_from_parallels(blade, trajectories[i])

        for j in range(0,len(trajectories[i])):
            list_joints =  map( lambda x: get_joint(base_num,joint_base,x), trajectories[i][j] )
            w_list, alpha_list, times = planning.compute_general_velocities(turb, list_joints, list_rays[j])
            torque = []
            sensible = []
            for ray, joints, w, alpha in izip(list_rays[j], list_joints, w_list, alpha_list):
                torque += [torque_computation(turbine, joints, w, alpha)]
                sensible += [sensibility(turbine, ray, w, alpha)]
                
            velocity_times += [(w_list, alpha_list, times)]
            torque_seg += [torque]
            sensible_seg += [sensible]
            
        velocities_times += [velocity_times]
        torque_traj += [torque_seg]
        sensible_traj += [sensible_seg]
    return velocities_times, torque_traj, sensible_traj

def check_angular_velocities(velocities):
    maxVel = turb.robot.GetDOFVelocityLimits()
    return all((abs(velocities)-maxVel)<0,1)

def check_angular_velocities_segs(velocity_segs, base_num):
    seg_base = DB.load_db_pickle(join(segs_path,str(base_num)+'.pkl'))
    joint_base = DB.load_db_pickle(join(joints_path,str(base_num)+'.pkl'))
    trajectories = seg_base[base_num]
    new_seg_base = dict()
    new_seg_base[base_num] = list()
    
    for i,velocity_seg in enumerate(velocity_segs):
        if len(velocity_seg)>0:
            new_seg = []
            for j,vel_alpha in enumerate(velocity_seg):
                vel = vel_alpha[0]
                check = check_angular_velocities(array(vel))
                indices = nonzero(check[1:] != check[:-1])[0] + 1
                b = split(array(trajectories[i][j]), indices)
                b = b[0::2] if check[0] else b[1::2]
                for bi in b:
                    if len(bi)>=3:
                        new_seg.append(list(bi))
            if len(new_seg)>0:
                new_seg_base[base_num].append(new_seg)
    return new_seg_base

def compute_new_segs(turb):
    new_path_seg = join(directory,'not_merged','new_seg')

    try:
        makedirs(new_path_seg)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
    
    while True:
        with open(join(directory,'fixed_db','db_visited_bases.pkl'), 'rb') as f:
            db_visited_bases = cPickle.load(f)
        base_num = None
        for key, value in db_visited_bases.iteritems():
            if value == False:
                db_visited_bases[key] = True
                base_num = key
                break
        if base_num is None:
            break
        with open(join(directory,'fixed_db','db_visited_bases.pkl'), 'wb') as f:
            cPickle.dump(db_visited_bases, f, cPickle.HIGHEST_PROTOCOL)
        del db_visited_bases

        try:
            seg_base = DB.load_db_pickle(join(segs_path,str(base_num)+'.pkl'))
        except IOError:
            continue

        velocities,_,_ = compute_velocities(base_num,turb)
        new_seg_base = check_angular_velocities_segs(velocities, base_num)

        print 'saving base_num: ', base_num

        try:
            DB.save_db_pickle(new_seg_base, join(new_path_seg,str(base_num)+'.pkl'))
        except IOError:
            raise 'Error saving db_base_to_joints.pkl'
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
    

if __name__ == '__main__':
    segs_path = 'db_lip/not_merged/seg'
    joints_path = 'db_lip/not_merged/joints'
    new_segs_path = 'db_lip/not_merged/new_seg'
    directory = 'db_lip'
    base_num = 13
    blade_folder = "lip"
    DB = db.DB(directory)

    dir_test = join(realpath('.'),'test')
    environ['OPENRAVE_DATA'] = str(dir_test)
    cfg = TurbineConfig.load('turbine_unittest.cfg','test')
    turb = Turbine(cfg)
    vis = Visualizer(turb.env)

    #compute_new_segs()
    velocities, torq, sens = compute_velocities(base_num,turb)
    #new_seg_base = check_angular_velocities_segs(velocities, base_num)
    #compute_new_segs(turb)
    #plot_segs_comp(base_num)
