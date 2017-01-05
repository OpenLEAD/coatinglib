#!/usr/bin/env python
import db
import blade_modeling
from turbine import Turbine
from turbine_config import TurbineConfig, ConfigFileError
import os
from visualizer import Visualizer
from os.path import join, isfile, realpath
from os import listdir
import rail_place
from numpy import array, random, zeros, dot, arange, linalg, sum
from numpy import sign, cross, random
from datetime import datetime
from os import makedirs
import cPickle
import errno
import mathtools
from math import pi
import mathtools
from copy import deepcopy
import ast

def merge():
    """
    Function to call merge method in db.py.
    It will merge the db files in 'not_merged' folder.
    """
    DB = db.DB(directory)
    DB.merge_db_directory(join(directory,'not_merged'))
    return

def plot_gradient():
    """
    Function to call plot_points_gradient in db.py.
    """
    DB = db.DB(directory)
    DB.plot_points_gradient(vis)
    DB.plot_bases_db(vis, turb)
    return

def generate_db():
    """
    Function to generate the db. It can be called multiple times
    by different process.
    """
    
    turb.robot.GetLink('Flame').Enable(False)
    blade = load_blade(blade_folder)
    DB = db.DB(directory)

    path = join(directory,'not_merged')
    try:
        makedirs(path)
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

        db_bases_to_num = DB.load_db_bases_to_num()
        for key, value in db_bases_to_num.iteritems():
            if value == base_num:
                base = key
                break
        del db_bases_to_num

        rp = rail_place.RailPlace(base)
        turb.place_rail(rp)
        turb.place_robot(rp)

        if turb.check_rail_collision():
            continue
        if turb.check_robotbase_collision():
            continue

        database = DB.generate_db(turb, blade, base_num)
        name = rp.getXYZ(turb.config)
        name = [round(name[0],3), round(name[1],3)]
        name = str(name)
        name = name.replace(', ','_')
        name = name.replace('[','')
        name = name.replace(']','')
        print 'saving base local (x,y): ', name
        DB.save_db_pickle(database, join(path,name+'.pkl'))
    return

def generate_robot_positions(number_of_positions=1000):
    """
    Function to generate random positions for the base of the robot.
    The positions will be saved and a db_visited_bases will be saved (False to everything).
    """
    rp = rail_place.rand_rail(turb.config, number_of_positions)
    keys = [ tuple(p.getPSAlpha()) for p in rp]
    
    psalpha_dict = dict(zip(keys, range(0,len(keys))))
    visited_bases = dict(zip(range(0,len(keys)), zeros(len(keys), dtype=bool)))

    path = join(directory,'fixed_db')
    
    try:
        makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
        
    with open(join(path,'db_bases_to_num.pkl'), 'wb') as f:
        cPickle.dump(psalpha_dict, f, cPickle.HIGHEST_PROTOCOL)
    with open(join(path,'db_visited_bases.pkl'), 'wb') as f:
        cPickle.dump(visited_bases, f, cPickle.HIGHEST_PROTOCOL)
    return       

def create_db_with_blade():
    """
    Funciton uses the blade points to create an empty db.
    """
    
    blade = load_blade(blade_folder)
    DB = db.DB(directory, blade)
    del blade
    return

def save_meridians(meridians):
    with open(join(directory,'fixed_db','meridians.pkl'), 'wb') as f:
        cPickle.dump(meridians, f, cPickle.HIGHEST_PROTOCOL)
    return

def save_parallels(parallels):
    with open(join(directory,'fixed_db','parallels.pkl'), 'wb') as f:
        cPickle.dump(parallels, f, cPickle.HIGHEST_PROTOCOL)
    return

def load_blade(folder):
    """
    Function to load blade model given folder.
    """
    xml_trajectories_path = join(folder,"trajectory/trajectory.xml")
    blade = blade_modeling.BladeModeling(turb, turb.blades[0])
    blade.load_trajectory(xml_trajectories_path)
    return blade
   
def make_grid():
    """
    See make_grid method in db.py.
    """
    
    DB = db.DB(directory)
    blade = load_blade(blade_folder_full)
    return DB.make_grid(blade)

def blade_borders(meridians):
    """
    After the grid creation,
    Find borders of the jiraublade and add it to meridians.
    """
    
    blade = load_blade(blade_folder_full)
    N = len(blade.trajectories)
    neg_border, pos_border = blade.find_borders(17, N-37)
    neg_1 = []; neg_2 = []; pos_1 = []; pos_2 = []
    for neg_b in neg_border:
        neg_1 += [neg_b[0][0:3]]
        neg_2 += [neg_b[-1][0:3]]
    for pos_b in pos_border:
        pos_1 += [pos_b[0][0:3]]
        pos_2 += [pos_b[-1][0:3]]   
    meridians = meridians[0:5]+[neg_1]+[neg_2]+meridians[5:10]+\
                [pos_1]+[pos_2]+meridians[10:]
    return meridians

def create_db_grid():
    """
    After the grid creation, borders and analysis if it makes sense,
    a db_grid must be created. See create_db_grid in db.py.
    """
    
    blade = load_blade(blade_folder)
    DB = db.DB(directory)
    DB.create_db_grid(blade)
    return

def grid_pick():
    """
    After the db_grid creation, some grids may not make sense.
    Iteractive remove those grids with this function.
    """
    
    blade = load_blade(blade_folder)
    DB = db.DB(directory)
    db_grid_to_mp = DB.load_db_grid_to_mp()
    db_grid_to_mp_copy = deepcopy(db_grid_to_mp)
    db_grid_to_bases = DB.load_db_grid_to_bases()
    db_grid_to_trajectories = DB.load_db_grid_to_trajectories()
    meridians = DB.load_grid_meridian()
    parallels = DB.load_grid_parallel()
    
    for key, value in db_grid_to_mp.iteritems():
        rays = DB.compute_rays_from_parallels(blade, db_grid_to_trajectories[key][0],
                                              db_grid_to_trajectories[key][1])
        s = vis.plot(meridians[value[0][0]],'meridian')
        s = vis.plot(meridians[value[0][1]],'meridian')
        s = vis.plot(parallels[value[1][0]],'parallel')
        s = vis.plot(parallels[value[1][1]],'parallel')
        s = vis.plot_lists(rays,'points', color=(1,0,0))
        print 'number of bases = ', len(db_grid_to_bases[key])
        x = raw_input('Remove base ? (y,n)')
        if x=='y':
            db_grid_to_mp_copy.pop(key, None)
            db_grid_to_bases.pop(key, None)
            db_grid_to_trajectories.pop(key, None)
            print "Grid removed"
        vis.remove_points('meridian')
        vis.remove_points('parallel')
        vis.remove_points('points')

        try:
            DB.save_db_pickle(db_grid_to_mp_copy, join(DB.path,'fixed_db','db_grid_to_mp.pkl'))
        except IOError: None

        try:
            DB.save_db_pickle(db_grid_to_bases, join(DB.path,'fixed_db','db_grid_to_bases.pkl'))
        except IOError: None

        try:
            DB.save_db_pickle(db_grid_to_trajectories, join(DB.path,'fixed_db','db_grid_to_trajectories.pkl'))
        except IOError: None

    return

def grid_add():
    """
    After db_grid creation, some grids may be added.
    """
    
    blade = load_blade(blade_folder)
    DB = db.DB(directory)
    
    db_grid_to_mp = DB.load_db_grid_to_mp()
    db_grid_to_bases = DB.load_db_grid_to_bases()
    db_grid_to_trajectories = DB.load_db_grid_to_trajectories()
    meridians = DB.load_grid_meridian()
    parallels = DB.load_grid_parallel()
    
    s = vis.plot_lists(meridians,'meridians')
    s = vis.plot(meridians[0],'meridians',color=(0,1,0))
    s = vis.plot(parallels[0],'parallels',color=(0,1,0))
    s = vis.plot_lists(parallels,'parallels')
    
    for key, value in db_grid_to_mp.iteritems():
        rays = DB.compute_rays_from_parallels(blade, db_grid_to_trajectories[key][0],
                                              db_grid_to_trajectories[key][1])
        s = vis.plot_lists(rays,'points', color=tuple(random.rand(3)))

    while True:
        db_grid_to_mp = DB.load_db_grid_to_mp()
        db_grid_to_bases = DB.load_db_grid_to_bases()
        db_grid_to_trajectories = DB.load_db_grid_to_trajectories()
        
        x = raw_input('Add: [(m1,m2),(p1,p2)]')
        grid = ast.literal_eval(x)
        key = max(db_grid_to_mp.keys())+1
        trajectories_in_grid, border = DB.get_points_in_grid(
            blade,[meridians[grid[0][0]],meridians[grid[0][1]]],
            [parallels[grid[1][0]],parallels[grid[1][1]]])
        rays = DB.compute_rays_from_parallels(blade, trajectories_in_grid,border)
        vis.plot_lists(rays,'rays')
        x = raw_input('Save ? (y,n)')
        if x == 'y':
            bases = DB.get_bases_trajectories(trajectories_in_grid)
            db_grid_to_mp[key] = grid
            db_grid_to_bases[key] = bases 
            db_grid_to_trajectories[key] = [trajectories_in_grid, border]
            try:
                DB.save_db_pickle(db_grid_to_mp, join(DB.path,'fixed_db','db_grid_to_mp.pkl'))
            except IOError: None
            try:
                DB.save_db_pickle(db_grid_to_bases, join(DB.path,'fixed_db','db_grid_to_bases.pkl'))
            except IOError: None
            try:
                DB.save_db_pickle(db_grid_to_trajectories, join(DB.path,'fixed_db','db_grid_to_trajectories.pkl'))
            except IOError: None

def compute_points_to_remove(grid_num):
    blade = load_blade(blade_folder)
    DB = db.DB(directory)
    db_grid_to_trajectories = DB.load_db_grid_to_trajectories()
    trajectories, borders = db_grid_to_trajectories[grid_num]
    rays = DB.compute_rays_from_parallels(blade, trajectories, borders)
    rays_traj = DB.compute_rays_from_parallels(blade, trajectories)

    def get_ray(point):
        model = blade.select_model(point)
        df = model.df(point)
        df = df/linalg.norm(df)
        return array(list(point)+list(df))
    
    points_to_remove = []
    border_to_remove = []
    for i in range(0,len(rays_traj)):
        
        if len(borders[i][0])>=0:
            if dot(get_ray(borders[i][0])[3:6],[0,0,1])>=0.45:
                border_to_remove.append(borders[i][0])
                borders[i][0] = []

        if len(borders[i][-1])>=0:
            if dot(get_ray(borders[i][-1])[3:6],[0,0,1])>=0.45:
                border_to_remove.append(borders[i][-1])
                borders[i][-1] = []

        if len(rays_traj[i])==0:
            continue
        for ray in rays_traj[i]:
            if dot(ray[3:6],[0,0,1])>=0.45:
                points_to_remove.append(tuple(ray[0:3]))
                
    vis.plot_lists(rays, color=(1,0,0))
    vis.plot(border_to_remove)
    vis.plot(points_to_remove)
    return borders, points_to_remove

def remove_points_from_db(grid_num, new_border, points_to_remove):
    blade = load_blade(blade_folder)
    DB = db.DB(directory)
    db_grid_to_trajectories = DB.load_db_grid_to_trajectories()
    trajectories, _ = db_grid_to_trajectories[grid_num]
    db_grid_to_trajectories[grid_num] = [trajectories, new_border]
    try:
        DB.save_db_pickle(db_grid_to_trajectories, join(DB.path,'fixed_db','db_grid_to_trajectories.pkl'))
    except IOError: None
    DB.remove_point(points_to_remove)
    return 

def grid_verticalization(grid_num):
    blade = load_blade(blade_folder)
    DB = db.DB(directory)
    db_grid_to_trajectories = DB.load_db_grid_to_trajectories()

    trajectories, borders = db_grid_to_trajectories[grid_num]
    rays = DB.compute_rays_from_parallels(blade, trajectories, borders)
    return mathtools.trajectory_verticalization(rays)

if __name__ == '__main__':

    directory = 'db'
    blade_folder = "jiraublade_hd_filtered"
    blade_folder_full = "jiraublade_hd"
    try:
        makedirs(directory)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
              
    dir_test = join(realpath('.'),'test')
    os.environ['OPENRAVE_DATA'] = str(dir_test)
    cfg = TurbineConfig.load('turbine_unittest.cfg','test')
    turb = Turbine(cfg)
    vis = Visualizer(turb.env)
    
    #generate_robot_positions()
    #create_db_with_blade()
    #generate_db()
    #merge()
    #plot_gradient()

    #meridians, parallels = make_grid()
    #meridians = blade_borders(meridians)
    #save_meridians(meridians)
    #save_parallels(parallels)
    #vis.plot_lists(meridians,'meridians')
    #vis.plot_lists(parallels,'parallels')

    #create_db_grid()
    #grid_pick()
    #grid_add()
    #feasible_bases = grid_validation(1)

    #borders, points_to_remove = compute_points_to_remove(24)
    #remove_points_from_db(24, borders, points_to_remove)



    
