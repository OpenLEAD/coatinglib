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
from numpy import sign, cross
from datetime import datetime
from os import makedirs
import cPickle
import errno
import mathtools
from math import pi
import mathtools

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
            

def trajectories_plot(trajectories):
    for traj in trajectories:
        vis.plot(traj, 'traj',(1,0,0))

def see_base_plot():
    path = 'db/converted'
    vis = Visualizer(turb.env)
    onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]

    for afile in onlyfiles:
        DB.plot_points_db(DB.load_db_pickle(join(path,afile)),vis)
        x = raw_input('wait')
        vis.remove_points('points_db')

def generate():
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

    
def make_grid(number_of_meridians = 12, number_of_meridians = 5 ):

    blade = load_blade(blade_folder)
    parallel = blade.trajectories[int(len(blade.trajectories)/2)]
    meridians = blade.draw_meridians_old(parallel, 1e-3, number_of_meridians)
    
    blade_full = load_blade(blade_folder_full)
    
    number_of_meridians = 5
    parallels = []
    step = len(blade_full.trajectories)/(number_of_meridians+1)
    for i in arange(step,len(blade_full.trajectories),step):
        parallels.append(blade_full.trajectories[i])
    
    return meridians, parallels
    
def get_points_in_grid(meridian1, meridian2, parallel1, parallel2):
    blade = load_blade(blade_folder)

    parallel_index_1 = int((blade.trajectory_iter_surface._Rn -
                            linalg.norm(parallel1[0][0:3]))/
                           blade.trajectory_iter_surface.coatingstep)
    parallel_index_2 = int((blade.trajectory_iter_surface._Rn -
                            linalg.norm(parallel2[0][0:3]))/
                           blade.trajectory_iter_surface.coatingstep)
    init = min(parallel_index_1,parallel_index_2)
    end = max(parallel_index_1,parallel_index_2)+1

    trajectories_in_grid = []

    for i in range(init,end):
        trajectory_in_grid = []
        parallel = array(blade.trajectories[i])
        blade.trajectory_iter_surface.find_iter(parallel[0])
        p1, sorted_parallel1 = closest_meridian_point(meridian1, parallel, blade)
        p2, sorted_parallel2 = closest_meridian_point(meridian2, parallel, blade)
        parallel1 = mathtools.direction_in_halfplane(parallel,p1[3:6])
        parallel2 = mathtools.direction_in_halfplane(parallel,p2[3:6])
        p1, sorted_parallel1 = closest_meridian_point(meridian1, parallel1, blade)
        p2, sorted_parallel2 = closest_meridian_point(meridian2, parallel2, blade)
        index_left = get_first_left_meridian_point_index(parallel, sorted_parallel1, p1)
        index_right = get_first_right_meridian_point_index(parallel, sorted_parallel2, p2)

        trajectory_in_grid += [p1]
        if abs(index_right - index_left)%(len(parallel)-1) == 1:
            pass
        elif index_left <= index_right:
            trajectory_in_grid += list(parallel[index_left:index_right+1])
        else:
            trajectory_in_grid += list(parallel[index_left:]) + list(parallel[:index_right+1])
            
        trajectory_in_grid += [p2]
    
        vis.plot(trajectory_in_grid,'traj_in_grid',color=(1,0,0))
        vis.plot(p2,'traj_in_grid',color=(0,1,0))
        vis.plot(p1,'traj_in_grid',color=(0,0,1))
    trajectories_in_grid.append(trajectory_in_grid)
    return trajectories_in_grid   

def closest_meridian_point(meridian, parallel, blade):
    min_dist = 100
    closest_meridian_point = []
    sorted_parallel = []
    for meridian_point in meridian:
        dist = sum((parallel[:,0:3]-meridian_point[0:3])*
                   (parallel[:,0:3]-meridian_point[0:3]),1)
        if min(dist) <= min_dist:
            closest_meridian_point = meridian_point
            min_dist = min(dist)
            sorted_parallel = [x for (y,x) in sorted(zip(dist,parallel))]
    model = blade.select_model(closest_meridian_point)
    return mathtools.curvepoint(model,blade.trajectory_iter_surface,closest_meridian_point[0:3]), sorted_parallel

def get_first_left_meridian_point_index(parallel, sorted_parallel, meridian_point):
    tan = cross(meridian_point[3:6],
                meridian_point[0:3]/linalg.norm(meridian_point[0:3]))
    for point in sorted_parallel:
        if sign(dot(tan,meridian_point[0:3]-point[0:3])) == 1:
            return parallel.tolist().index(list(point))

def get_first_right_meridian_point_index(parallel, sorted_parallel, meridian_point):
    tan = cross(meridian_point[3:6],
                meridian_point[0:3]/linalg.norm(meridian_point[0:3]))

    for point in sorted_parallel:
        if sign(dot(tan,meridian_point[0:3]-point[0:3])) == -1:
            return parallel.tolist().index(list(point))

def compute_bases_to_coat_points(trajectories_in_grid):
    DB = db.DB(directory)
    bases = DB.compute_bases_to_coat_points(trajectories_in_grid)
    print len(bases)
    for base in bases:
        vis.plot(rail_place.RailPlace(base).getXYZ(turb.config), 'bases',(0,1,0))
    return bases

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
    
    #generate_robot_positions()
    #create_db_with_blade()

    #rename_files_by_base()
    #convert()
    #vis = Visualizer(turb.env)
    #generate()
    #merge()

    vis = Visualizer(turb.env)
    #plot_gradient()
    #blade_segmentation()
    #plot_points_covered_by_n(600,50)

    meridians, parallels = make_grid()
##    trajectories_in_grid = get_points_in_grid(meridians[6], meridians[7], parallels[4], parallels[5])
##    compute_bases_to_coat_points(trajectories_in_grid)
