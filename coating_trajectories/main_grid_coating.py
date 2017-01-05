from math import pi
import db
from os.path import join, realpath
import os
from turbine import Turbine
from turbine_config import TurbineConfig, ConfigFileError
from visualizer import Visualizer
from numpy import ones, array, dot, linspace
from openravepy import matrixFromAxisAngle
import rail_place
import mathtools
import blade_modeling
import planning
from copy import deepcopy
import csv

def base_for_grid_coating(grid_num):
    
    bases = []
    scores = []
    db_grid_to_trajectories = dict_angle_db[0].load_db_grid_to_trajectories()
    trajectories, borders = db_grid_to_trajectories[grid_num]
    for angle, DB in dict_angle_db.iteritems():
        base, score = DB.get_best_bases_trajectories(trajectories)
        bases.append(base)
        scores.append(score)

    angles_1 = []
    bases_1 = []
    scores_1 = []
    for i in range(0,len(scores)):
        angles_1 +=list(ones(len(scores[i]))*db_angles[i])
        bases_1 += list(bases[i])
        scores_1 += list(scores[i])
    sorted_bases = sorted(zip(-array(scores_1),bases_1,angles_1))
    
    return sorted_bases, trajectories, borders

def plot_base_for_grid_coating(sorted_bases, trajectories, borders):

    counter = 0
    while counter<len(sorted_bases):
        score, base, angle = sorted_bases[counter]
        
        blade = dict_angle_blade[angle]
        rborders = mathtools.rotate_points(borders, matrixFromAxisAngle([angle,0,0]))
        rays = dict_angle_db[angle].compute_rays_from_parallels(blade, trajectories, rborders)
        vis.plot_lists(rays,'rays')
        
        for blade in turb.blades:
            blade.SetTransform(matrixFromAxisAngle([angle,0,0]))

        rp = rail_place.RailPlace(base)
        turb.place_rail(rp)
        turb.place_robot(rp)
        print 'score: ', -score, '\t rotor angle: ', angle, '\t robot position:', rp.getXYZ(cfg)

        x = raw_input('next ? (y,n)')

        for blade in turb.blades:
            blade.SetTransform(matrixFromAxisAngle([-angle,0,0]))

        if x=='n':
            return score, base, angle
        vis.remove_points('rays')

        counter+=1
        
def load_blade(folder):
    """
    Function to load blade model given folder.
    """
    xml_trajectories_path = join(folder,"trajectory/trajectory.xml")
    blade = blade_modeling.BladeModeling(turb, turb.blades[0])
    blade.load_trajectory(xml_trajectories_path)
    return blade       

def base_grid_validation(blade_angle, rays):
    """
    Given the real blade angle:
    1) rotate the blades (update the environment);
    2) rotate the RBF model of the blade;
    3) rotate grid points;
    4) organize trajectories (removing empty points, adding borders,
    and making one full zigzagging list);
    5) compute optimization.

    Keyword arguments:
    blade_angle -- real angle of the blade
    rays -- points to be coated (list-n of lists-mx6) w.r.t. the 0 rotor angle.
    """
    
    T = matrixFromAxisAngle([blade_angle,0,0])
    for blade in turb.blades:
        blade.SetTransform(T)
        
    blade = dict_angle_blade[0]
    blade.rotate_models(T)

    rotated_rays = deepcopy(rays)

    rotated_rays = mathtools.rotate_trajectories(rotated_rays, T)
    organized_rays = []
    for i in range(0,len(rotated_rays)):
        if i%2==0:
            organized_rays += rotated_rays[i]
        else:
            organized_rays += reversed(rotated_rays[i])
    organized_rays = [x for x in organized_rays if x != []]

    turb.robot.GetLink('Flame').Enable(False)

    joint_solutions = planning.compute_robot_joints_opt(turb, organized_rays, 0,
                                                    blade.trajectory_iter_surface)
    score = len(joint_solutions)*1.0/len(organized_rays)
    
    return joint_solutions, score

def tolerance_test(sorted_bases, trajectories, borders, score_threshold=0.9):
    
    blade = dict_angle_blade[0]
    rays = dict_angle_db[0].compute_rays_from_parallels(blade, trajectories, borders)

    table = ['Rotor Angle', 'Base Position', 'Expected Score', 'Score']
    tolerance_list = []

    try:
        with open(join('tolerance_db',str(grid_num)+'_db.csv'),'rb') as f:
            reader = csv.reader(f)
            next(reader,None)
            tolerance_list = list(reader)
    except IOError:
        with open(join('tolerance_db',str(grid_num)+'_db.csv'),'w') as f:
            writer = csv.writer(f)
            writer.writerow(table)

    for i in range(0,len(tolerance_list)):
        tolerance_list[i] = map(float,tolerance_list[i])
        tolerance_list[i][0:2] = map(int,tolerance_list[i][0:2])
    
    counter = 0
    for i in range(0,len(sorted_bases)):

        expected_score, base, angle = sorted_bases[i]
        if -expected_score < score_threshold: break

        rp = rail_place.RailPlace(base)
        turb.place_rail(rp)
        turb.place_robot(rp)
        base_num = (db.DB(db_directories[0])).load_db_bases_to_num()[base]

        if len(tolerance_list)>0:
            if base_num in (array(tolerance_list)[:,1]).tolist():
                continue

        for ang in linspace(angle - 10*pi/180, angle + 10*pi/180, 21):
            joint_solutions, score = base_grid_validation(ang, rays)
            tolerance_list.append([round(ang*180/pi), base_num, -expected_score, score])

        print('iter %3i ' % (counter))
        counter+=1
        with open(join('tolerance_db',str(grid_num)+'_db.csv'),'w') as f:
            writer = csv.writer(f)
            writer.writerow(table)
            writer.writerows(tolerance_list)
    
    return

def non_coatable_grids():
    non_coatable = []
    for grid, trajs in dict_angle_db[0].load_db_grid_to_trajectories().iteritems():
        sorted_bases, trajectories, borders = base_for_grid_coating(grid)
        expected_score, base, angle = sorted_bases[0]
        if -expected_score != 1:
            non_coatable.append(grid)
    return non_coatable

if __name__ == '__main__':

    #-----------------------------------------------------------------------
    # DB inputs
    db_directories = ['db']#, 'db_45', 'db_-45']
    db_angles = [0]#, pi/4, -pi/4]
    blade_folder = ['jiraublade_hd_filtered']#, 'jiraublade_hd_filtered_45',
                   # 'jiraublade_hd_filtered_-45']
    #-----------------------------------------------------------------------
    
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
    
    #vis = Visualizer(turb.env)

    # Grid input
    grid_num = 75
    sorted_bases, trajectories, borders = base_for_grid_coating(grid_num)
    #plot_base_for_grid_coating(sorted_bases, trajectories, borders)
    #tolerance_test(sorted_bases, trajectories, borders, init=63)
