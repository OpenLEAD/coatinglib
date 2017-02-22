#!/usr/bin/env python
import db
from os.path import join, realpath
from os import environ
import blade_modeling
from turbine import Turbine
from turbine_config import TurbineConfig, ConfigFileError
import planning
from visualizer import Visualizer

def load_blade(folder):
    """
    Function to load blade model given folder.
    """
    xml_trajectories_path = join(folder,"trajectory/trajectory.xml")
    blade = blade_modeling.BladeModeling(turb, turb.blades[0])
    blade.load_trajectory(xml_trajectories_path)
    return blade

def compute_velocities():
    DB = db.DB(directory)
    blade = load_blade(blade_folder)
    
    seg_base = DB.load_db_pickle(join(segs_path,str(base_num)+'.pkl'))
    joint_base = DB.load_db_pickle(join(joints_path,str(base_num)+'.pkl'))
    trajectories = seg_base[base_num]
    velocities_times = []

    def get_joint(point_num):
        return joint_base[base_num][point_num]
    
    for i in range(0,len(trajectories)):
        velocity_times = []
        if len(trajectories[i])>0:
            list_rays = DB.compute_rays_from_parallels(blade, trajectories[i])

        for j in range(0,len(trajectories[i])):
            list_joints =  map( lambda x: get_joint(x), trajectories[i][j] ) 
            velocity_times.append(planning.compute_general_velocities(turb, list_joints, list_rays[j]))
        velocities_times.append(velocity_times)
    return velocities_times

def check_angular_velocities(velocities):
    maxVel = turb.robot.GetDOFMaxVel()
    for velocity in velocities:
        if any((abs(velocity)-maxVel)>0):
            return False
    return True

def check_angular_velocities_segs(velocity_segs):
    check_lists = []
    for velocity_seg in velocity_segs:
        check_list = []
        if len(velocity_seg)>0:
            check_list = []
            for vel_alpha in velocity_seg:
                vel = vel_alpha[0]
                check_list.append(check_angular_velocities(vel))
        check_lists.append(check_list)
    return check_lists

if __name__ == '__main__':
    segs_path = 'db/not_merged/seg'
    joints_path = 'db/not_merged/joints'
    directory = 'db'
    base_num = 0
    blade_folder = "jiraublade_hd_filtered"

    dir_test = join(realpath('.'),'test')
    environ['OPENRAVE_DATA'] = str(dir_test)
    cfg = TurbineConfig.load('turbine_unittest.cfg','test')
    turb = Turbine(cfg)
    vis = Visualizer(turb.env)

    velocities = compute_velocities()
    check_lists = check_angular_velocities_segs(velocities)
