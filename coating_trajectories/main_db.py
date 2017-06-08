#!/usr/bin/env python
import db
import blade_modeling
from turbine import Turbine
from turbine_config import TurbineConfig, ConfigFileError
import os
from visualizer import Visualizer
from os.path import join, isfile, realpath, basename
from os import listdir
import rail_place
from numpy import array, random, zeros, dot, arange, linalg, sum
from numpy import sign, cross, random
from datetime import datetime
from os import makedirs
import cPickle
import errno
from math import pi, tan, atan2 
import mathtools
from copy import deepcopy
import ast
from openravepy import matrixFromAxisAngle
import sensibility
import blade_coverage

def generate_db_joints():
    """
    Function to generate the trajectory_db and joints_db.
    It can be called multiple times by different process.
    """
    
    turb.robot.GetLink('Flame').Enable(False)
    path_joints = join(path,'joints')
    path_seg = join(path,'seg')

    try:
        makedirs(path_joints)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

    try:
        makedirs(path_seg)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
    
    while True:
        visited_bases = db.load_pickle(join(path,'visited_bases.pkl'))
        base_num = None
        for key, value in visited_bases.iteritems():
            if value == False:
                visited_bases[key] = True
                base_num = key
                break
        if base_num is None:
            break
        db.save_pickle(visited_bases,join(path,'visited_bases.pkl'))
        del visited_bases

        base_to_joints, base_to_seg = DB.generate_db_joints(
            base_num, minimal_number_of_points_per_trajectory, do_side_filter)
        if base_to_joints == None:
            continue
        
        print 'saving base_num: ', base_num

        try:
            db.save_pickle(base_to_joints, join(path_joints,str(base_num)+'.pkl'))
        except IOError:
            raise 'Error saving db_base_to_joints.pkl'
        
        try:
            db.save_pickle(base_to_seg, join(path_seg,str(base_num)+'.pkl'))
        except IOError:
            raise 'Error saving db_base_to_seg.pkl'
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

    db.save_pickle(psalpha_dict,join(path,'bases_to_num.pkl'))
    db.save_pickle(visited_bases,join(path,'visited_bases.pkl'))
    return       

def create_db():
    """
    Funciton uses the blade points to create an empty db.
    """
    
    DB.create_db()
    return

def make_grid(number_of_meridians, number_of_parallels, init_parallel):
    """
    See make_grid method in db.py.
    """   
    return DB.make_grid(number_of_meridians, number_of_parallels, init_parallel)

def blade_borders(meridians):
    """
    After the grid creation,
    Find borders of the jiraublade and add it to meridians.
    """
    
    blade = DB.load_blade_full()
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
    DB.create_db_grid()
    return

def grid_pick(vis):
    """
    After the db_grid creation, some grids may not make sense.
    Iteractive remove those grids with this function.
    """

    grid_to_mp = DB.load_grid_to_mp()
    grid_to_mp_copy = deepcopy(grid_to_mp)
    grid_to_trajectories = DB.load_grid_to_trajectories()
    meridians = DB.load_grid_meridian()
    parallels = DB.load_grid_parallel()
    
    for key, value in grid_to_mp.iteritems():
        rays = DB.compute_rays_from_parallels(grid_to_trajectories[key][0],
                                              grid_to_trajectories[key][1])
        s = vis.plot(meridians[value[0][0]],'meridian')
        s = vis.plot(meridians[value[0][1]],'meridian')
        s = vis.plot(parallels[value[1][0]],'parallel')
        s = vis.plot(parallels[value[1][1]],'parallel')
        s = vis.plot_lists(rays,'points', color=(1,0,0))
        x = raw_input('Remove base ? (y,n)')
        if x=='y':
            grid_to_mp_copy.pop(key, None)
            grid_to_trajectories.pop(key, None)
            print "Grid removed"
        vis.remove_points('meridian')
        vis.remove_points('parallel')
        vis.remove_points('points')

        try:
            db.save_pickle(grid_to_mp_copy, join(DB.path,'grid_to_mp.pkl'))
        except IOError: None

        try:
            db.save_pickle(grid_to_trajectories, join(DB.path,'grid_to_trajectories.pkl'))
        except IOError: None

    return

def grid_add(vis):
    """
    After db_grid creation, some grids may be added.
    """
    
    grid_to_mp = DB.load_grid_to_mp()
    grid_to_trajectories = DB.load_grid_to_trajectories()
    meridians = DB.load_grid_meridian()
    parallels = DB.load_grid_parallel()
    
    s = vis.plot_lists(meridians,'meridians')
    s = vis.plot(meridians[0],'meridians',color=(0,1,0))
    s = vis.plot(parallels[0],'parallels',color=(0,1,0))
    s = vis.plot_lists(parallels,'parallels')
    
    for key, value in grid_to_mp.iteritems():
        rays = DB.compute_rays_from_parallels(grid_to_trajectories[key][0],
                                              grid_to_trajectories[key][1])
        s = vis.plot_lists(rays,'points', color=tuple(random.rand(3)))

    while True:
        grid_to_mp = DB.load_grid_to_mp()
        grid_to_trajectories = DB.load_grid_to_trajectories()
        
        x = raw_input('Add: [(m1,m2),(p1,p2)]')
        grid = ast.literal_eval(x)
        key = max(grid_to_mp.keys())+1
        trajectories_in_grid, border = DB.get_points_in_grid(
            [meridians[grid[0][0]],meridians[grid[0][1]]],
            [parallels[grid[1][0]],parallels[grid[1][1]]])
        rays = DB.compute_rays_from_parallels(trajectories_in_grid,border)
        vis.plot_lists(rays,'rays')
        x = raw_input('Save ? (y,n)')
        if x == 'y':
            bases = DB.get_bases_trajectories(trajectories_in_grid)
            grid_to_mp[key] = grid
            grid_to_trajectories[key] = [trajectories_in_grid, border]
            try:
                db.save_pickle(grid_to_mp, join(DB.path,'grid_to_mp.pkl'))
            except IOError: None
            try:
                db.save_pickle(grid_to_trajectories, join(DB.path,'grid_to_trajectories.pkl'))
            except IOError: None

def compute_points_to_remove(grid_num, vis):

    blade = DB.load_blade()
    grid_to_trajectories = DB.load_grid_to_trajectories()
    trajectories, borders = grid_to_trajectories[grid_num]
    rays = DB.compute_rays_from_parallels(trajectories, borders)
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

    grid_to_trajectories = DB.load_grid_to_trajectories()
    trajectories, _ = grid_to_trajectories[grid_num]
    grid_to_trajectories[grid_num] = [trajectories, new_border]
    try:
        db.save_pickle(grid_to_trajectories, join(DB.path,'grid_to_trajectories.pkl'))
    except IOError: None
    DB.remove_point(points_to_remove)
    return 

def clear_visited_bases():
    DB.clear_visited_bases()
    return

def create_db_from_segments(directory):
    db_main = DB.create_db_from_segments(directory)
    db.save_pickle(db_main, join(path,'db.pkl') )
    return 

def turbine_rotate(turb):
    for blade in turb.blades:
        T = blade.GetTransform()
        blade.SetTransform(dot(DB.T,T))
    return

def bases_grids_coating(threshold):
    
    grid_to_trajectories = DB.load_grid_to_trajectories()
    grid_bases = dict()
    for grid_num in grid_to_trajectories.keys():
        bases, scores = DB.base_grid_coating(grid_num)
        xy = []
        for i, base in enumerate(bases):
            if scores[i] >= threshold:
                rp = rail_place.RailPlace(base)
                xyz = rp.getXYZ(cfg)
                value = set([(xyz[0],xyz[1])])
                grid_bases[grid_num] = grid_bases.get(grid_num,set()) | value
            else: break
    db.save_pickle(grid_bases,join(DB.db_main_path,'grid_bases.pkl'))
    return

def verify_base_grid_threshold():
    grid_to_trajectories = DB.load_grid_to_trajectories()
    grid_bases = dict()
    for grid_num in grid_to_trajectories.keys():
        trajectories, borders = grid_to_trajectories[grid_num]
        bases, scores = DB.base_grid_coating(grid_num)
        print 'Grid: ', grid_num, '/ score: ', scores[0]
    return

def plot_grid_coat(vis):
    bases, scores = DB.base_grid_coating(grid_num)
    non_coatable = DB.plot_grid_coat(vis, grid_num,
                                     DB.load_bases_to_num()[bases[0]])
    return

def _compute_lines(x_step=.05, angle_step=pi/18):
    lines=[]
    x_range = [turb.config.environment.x_min,turb.config.environment.x_max]
    angle_range = [-pi/2,pi/2]
    for x in arange(x_range[0],x_range[1],x_step):
        for angle in arange(angle_range[0]+1,angle_range[1],angle_step):
            tanalpha = tan(angle)
            lines.append(((x,0),(x-tanalpha,1.)))
    lines.append(((0,0),(1,0)))
    return lines

def generate_rail_configurations():
    threshold = 0.2
    lines = _compute_lines()
    DB.compute_rail_configurations(lines, threshold)
    return

def make_validate_file():
    
    grids = []
    gr = DB.info.findall('grids')
    for g in gr:
        grids+=eval('blade_coverage.'+g.text+'()')
    grids = list(set(grids))

    db_grids = DB.get_dbs_grids()
    grid_db_lines = dict()
    for grid in grids:
        db_lines = dict()
        for dbi in db_grids.keys():
            if grid in db_grids[dbi]:
                line_grid = db.load_pickle(join(dbi,'rails','line_grid.pkl'))
                line_grid_dist = db.load_pickle(join(dbi,'rails','line_grid_dist.pkl'))
                visited_lines = dict()
                lines = line_grid.keys()
                for line in lines:
                    if grid in line_grid[line]:
                        visited_lines[line] = False
                db_lines[dbi] = visited_lines
        grid_db_lines[grid] = db_lines
    db.save_pickle(grid_db_lines,join(DB.path,'visited_lines.pkl'))

def validate_bases():
    """
    Function to validate the bases generated.
    It can be called multiple times by different process.
    """
    
    turb.robot.GetLink('Flame').Enable(False)
    
    while True:
        visited_lines = db.load_pickle(join(DB.path,'visited_lines.pkl'))
        grids = visited_lines.keys()
        for grid in grids:
            dbs = visited_lines[grid].keys()
            for dbi in dbs:
                lines = visited_lines[grid][dbi].keys()
                for line in lines:
                    if visited_lines[grid][dbi][line] == False:
                        visited_lines[grid][dbi][line] = True
                        db_to_test = dbi
                        line_to_test = line
                        grid_to_test = grid
                        break
                else: continue
                break
            else: continue
            break
        else: break
        
        db.save_pickle(visited_lines,join(path,'visited_lines.pkl'))
        del visited_lines

        dbs = DB.info.findall('db')
        for dbi in dbs:
            if dbi.find('name').text == basename(db_to_test):
                DB.T = DB._extract_T(dbi)
                DB.db_main_path = join(DB.path,dbi.find('path').text)
        line_grid = db.load_pickle(join(db_to_test,'rails','line_grid.pkl'))
        line_grid_dist = db.load_pickle(join(db_to_test,'rails','line_grid_dist.pkl'))
        try:
            point_near, distance, distance_str = line_grid_dist[line_to_test][grid_to_test]
        except KeyError:
            continue
        x1 = line_to_test[0][0]; y1 = line_to_test[0][1]
        x2 = line_to_test[1][0]; y2 = line_to_test[1][1]
        p = mathtools.closest_point_line_3d(array(line_to_test[0]), array(line_to_test[1]), point_near)
        psa = (x1, sign(p[1])*linalg.norm(p-line_to_test[0]), sign(p[1])*atan2(x1-p[0],abs(p[1]-y1)))
        res = blade_coverage.base_grid_validation(turb, psa, DB, grid)
        if not res.success:
            line_grid[line_to_test] = line_grid[line_to_test] - set([grid])
            db.save_pickle(line_grid,join(db_to_test,'rails','line_grid.pkl'))
            _ = line_grid_dist[line_to_test].pop(grid,None)
            db.save_pickle(line_grid_dist,join(db_to_test,'rails','line_grid_dist.pkl'))
    return

if __name__ == '__main__':

    area = 'FACE'
    db_name = ''
    path = join(area,db_name)
    
    dir_test = join(realpath('.'),'test')
    os.environ['OPENRAVE_DATA'] = str(dir_test)
    cfg = TurbineConfig.load('turbine_unittest.cfg','test')
    turb = Turbine(cfg)
    
    DB = db.DB(area, turb, db_name)
    turbine_rotate(turb)
    
    threshold = 0.90
    grid_num = 2
    
    #vis = Visualizer(turb.env)
    
    #generate_robot_positions()
    #create_db()
    #clear_visited_bases()

    # Generate robot joints inputs
    minimal_number_of_points_per_trajectory, do_side_filter = 3, False
    #generate_db_joints()

    #meridians, parallels = make_grid(number_of_meridians = 4, number_of_parallels = 2, init_parallel = 0)
    #meridians = blade_borders(meridians)

    #create_db_grid()
    #grid_pick()
    #grid_add()

    #borders, points_to_remove = compute_points_to_remove(24)
    #remove_points_from_db(24, borders, points_to_remove)

    #sensibility.compute_new_segs_DB(DB)
    #create_db_from_segments(join(path,'new_seg'))

    #plot_grid_coat()
    #verify_base_grid_threshold()
    #bases_grids_coating(threshold)
    #DB.plot_convex_grid(threshold,grid_num)

    #generate_rail_configurations()

    #make_validate_file()
    #validate_bases()
    
