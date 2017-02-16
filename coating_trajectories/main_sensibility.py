#!/usr/bin/env python
import db
from os.path import join, realpath
from os import environ
import blade_modeling
from turbine import Turbine
from turbine_config import TurbineConfig, ConfigFileError
import planning

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
    velocities = []

    def get_joint(point_num):
        return joint_base[base_num][point_num]
    
    for i in range(0,len(trajectories)):
        velocity = []
        if len(trajectories[i])>0:
            list_rays = DB.compute_rays_from_parallels(blade, trajectories[i])

        for j in range(0,len(trajectories[i])):
            list_joints =  map( lambda x: get_joint(x), trajectories[i][j] ) 
            velocity.append(planning.compute_general_velocities(turb, list_joints, list_rays[j]))
        velocities.append(velocity)
    return velocities
                         

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

    velocities = compute_velocities()
