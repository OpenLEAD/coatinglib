from math import pi
import db
from os.path import join, realpath
import os
from turbine import Turbine
from turbine_config import TurbineConfig, ConfigFileError
from visualizer import Visualizer
from numpy import ones, array, dot
from openravepy import matrixFromAxisAngle
import rail_place
import mathtools
import blade_modeling

def base_for_grid_coating(grid_num):
    
    bases = []
    scores = []
    db_grid_to_trajectories = dict_angle_db[0].load_db_grid_to_trajectories()
    trajectories, borders = db_grid_to_trajectories[grid_num]
    for key, value in dict_angle_db.iteritems():
        base, score = value.get_best_bases_trajectories(trajectories)
        bases.append(base)
        scores.append(score)
    return bases, scores

def plot_base_for_grid_coating(grid_num, bases, scores):

    db_grid_to_trajectories = dict_angle_db[0].load_db_grid_to_trajectories()
    trajectories, borders = db_grid_to_trajectories[grid_num]
    
    angles_1 = []
    bases_1 = []
    scores_1 = []
    for i in range(0,len(scores)):
        angles_1 +=list(ones(len(scores[i]))*db_angles[i])
        bases_1 += list(bases[i])
        scores_1 += list(scores[i])
    sorted_bases = sorted(zip(-array(scores_1),bases_1,angles_1))

    counter = 0
    while counter<len(sorted_bases):
        score, base, angle = sorted_bases[counter]
        print 'score: ', -score, '\t angle: ', angle
        
        blade = dict_angle_blade[angle]
        rborders = mathtools.rotate_points(borders, matrixFromAxisAngle([angle,0,0]))
        rays = dict_angle_db[angle].compute_rays_from_parallels(blade, trajectories, rborders)
        vis.plot_lists(rays,'rays')
        
        for blade in turb.blades:
            blade.SetTransform(matrixFromAxisAngle([angle,0,0]))

        rp = rail_place.RailPlace(base)
        turb.place_rail(rp)
        turb.place_robot(rp)

        x = raw_input('next ? (y,n)')
        if x=='n':
            return score, base, angle

        vis.remove_points('rays')
        for blade in turb.blades:
            blade.SetTransform(matrixFromAxisAngle([-angle,0,0]))

        counter+=1
        
def load_blade(folder):
    """
    Function to load blade model given folder.
    """
    xml_trajectories_path = join(folder,"trajectory/trajectory.xml")
    blade = blade_modeling.BladeModeling(turb, turb.blades[0])
    blade.load_trajectory(xml_trajectories_path)
    return blade       

if __name__ == '__main__':

    db_directories = ['db', 'db_45', 'db_-45']
    db_angles = [0, pi/4, -pi/4]
    blade_folder = ['jiraublade_hd_filtered', 'jiraublade_hd_filtered_45',
                    'jiraublade_hd_filtered_-45']
    
    dir_test = join(realpath('.'),'test')
    os.environ['OPENRAVE_DATA'] = str(dir_test)
    cfg = TurbineConfig.load('turbine_unittest.cfg','test')
    turb = Turbine(cfg)

    dict_angle_db = dict()
    dict_angle_blade = dict()
    for i in range(0,len(db_directories)):
        dict_angle_db[db_angles[i]]=db.DB(db_directories[i])
    for i in range(0,len(blade_folder)):
        dict_angle_blade[db_angles[i]]=load_blade(blade_folder[i])
    
    vis = Visualizer(turb.env)

    grid_num = 0

    bases, scores = base_for_grid_coating(grid_num)
    plot_base_for_grid_coating(grid_num, bases, scores)
