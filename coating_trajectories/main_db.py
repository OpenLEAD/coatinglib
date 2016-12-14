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
    DB = db.DB(directory)
    DB.merge_db_directory(join(directory,'not_merged'))
    return

def plot_gradient():
    DB = db.DB(directory)
    DB.plot_points_gradient(vis)
    DB.plot_bases_db(vis, turb)
    return

def plot_points_covered_by_n(n_max, n_min=0):
    DB = db.DB(directory)
    points = []
    data = DB.load_db()
    db_points_to_num = DB.load_db_points_to_num()
    all_points_tuple, all_points_num = db_points_to_num.keys(), db_points_to_num.values()
    del db_points_to_num
    all_points_tuple = [x for (y,x) in sorted(zip(all_points_num,all_points_tuple))]
    del all_points_num
    for key, value in data.iteritems():
        if len(value) <= n_max and len(value) >= n_min:
            points.append(all_points_tuple[key])
    vis.plot(points,'plot_points_covered_by_n')
    return

def generate_db():
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

def generate_robot_positions():
    rp = rail_place.rand_rail(turb.config, 1000)
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
    blade = load_blade(blade_folder)
    DB = db.DB(directory, blade)
    del blade
    return

def load_meridians():
    with open(join(directory,'fixed_db','meridians.pkl'), 'rb') as f:
            return cPickle.load(f)

def save_meridians(meridians):
    with open(join(directory,'fixed_db','meridians.pkl'), 'wb') as f:
        cPickle.dump(meridians, f, cPickle.HIGHEST_PROTOCOL)
    return

def load_parallels():
    with open(join(directory,'fixed_db','parallels.pkl'), 'rb') as f:
            return cPickle.load(f)

def save_parallels(parallels):
    with open(join(directory,'fixed_db','parallels.pkl'), 'wb') as f:
        cPickle.dump(parallels, f, cPickle.HIGHEST_PROTOCOL)
    return

def load_blade(folder):
    xml_trajectories_path = join(folder,"trajectory/trajectory.xml")
    blade = blade_modeling.BladeModeling(turb, turb.blades[0])
    blade.load_trajectory(xml_trajectories_path)
    return blade
   
def make_grid():
    DB = db.DB(directory)
    blade = load_blade(blade_folder_full)
    return DB.make_grid(blade)

def compute_bases_to_coat_points(trajectories_in_grid):
    DB = db.DB(directory)
    bases = DB.get_bases_trajectories(trajectories_in_grid)
    print len(bases)
    for base in bases:
        vis.plot(rail_place.RailPlace(base).getXYZ(turb.config), 'bases',(0,1,0))
    return bases

def blade_borders(meridians):
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

def get_points_in_grid(meridian1, meridian2, parallel1, parallel2):
    blade = load_blade(blade_folder)
    DB = db.DB(directory)
    return DB.get_points_in_grid(blade, [meridian1, meridian2], [parallel1, parallel2])

def create_db_grid():
    blade = load_blade(blade_folder)
    meridians, parallels = load_meridians(), load_parallels()
    DB = db.DB(directory)
    DB.create_db_grid(blade, meridians, parallels)
    return
   
if __name__ == '__main__':

    directory = 'new_db'
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
    #meridians, parallels = load_meridians(), load_parallels()
    #meridians = blade_borders(meridians)
    #save_meridians(meridians)
    #save_parallels(parallels)
    #vis.plot_lists(meridians,'meridians')
    #vis.plot_lists(parallels,'parallels')

    db_grid_to_mp, db_grid_to_bases, db_grid_to_trajectories = create_db_grid()


    
