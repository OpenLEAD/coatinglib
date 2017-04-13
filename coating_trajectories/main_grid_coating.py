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
from copy import deepcopy, copy
import csv
import cPickle 

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

def grid_verticalization(grid_num):
    DB = dict_angle_db[db_angles[0]]
    blade = dict_angle_blade[db_angles[0]]
    db_grid_to_trajectories = DB.load_db_grid_to_trajectories()

    trajectories, borders = db_grid_to_trajectories[grid_num]
    rays = DB.compute_rays_from_parallels(blade, trajectories, borders)
    return mathtools.trajectory_verticalization(rays)

def compute_best_intersection(intersection, grid_nums, grid_bases):

    best_intersection_len = 0
    best_intersection = None
    intersections = None
    for grid in grid_nums:
        intersections = mathtools.intersection(intersection, grid_bases[grid], 0.2)
        if len(intersections)>best_intersection_len:
            best_intersection_len = len(intersections)
            best_intersection = grid
            intersection = intersections
    return best_intersection, best_intersection_len, intersection

def compute_best_intersections(grid_num, grid_nums, min_intersections=2, grid_bases=None):

    if grid_bases==None:
        grid_bases = bases_for_all_grids(threshold)

    all_grids = copy(grid_nums)
    remove_grids = [grid_num]
    intersection = grid_bases[grid_num]
    while True:
        for grid in remove_grids:
            try:
               grid_nums.remove(grid)
            except: None
            
        best_intersection, best_intersection_len, intersection = compute_best_intersection(
            intersection, grid_nums, grid_bases)

        if best_intersection_len<=min_intersections:
            return intersection, remove_grids

        remove_grids.append(best_intersection)
        grid_nums = copy(all_grids)

        if len(remove_grids) == len(all_grids):
            break
    return intersection, remove_grids
