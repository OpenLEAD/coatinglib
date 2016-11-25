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
from numpy import array, random, zeros
from datetime import datetime
from os import makedirs
import cPickle
import errno

def convert():
    DB.convert_db_point_base_directory('db/servidor', 'db/converted')

def merge():
    DB = db.DB(directory)
    DB.merge_db_directory(join(directory,'not_merged'))

def plot_gradient():
    DB = db.DB(directory)
    DB.plot_points_gradient(vis)
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
    xml_trajectories_path = join(blade_folder,"trajectory/trajectory.xml")
    blade = blade_modeling.BladeModeling(turb, turb.blades[0])
    blade.load_trajectory(xml_trajectories_path)
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
        name = rp.getXYZ(turb)
        name = [round(name[0],3), round(name[1],3)]
        name = str(name)
        name = name.replace(', ','_')
        name = name.replace('[','')
        name = name.replace(']','')
        print 'saving base local (x,y): ', name
        DB.save_db_pickle(database, join(path,name+'.pkl'))
    return

def rename_files_by_base():
    path = 'db/converted'
    onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]

    db_bases_to_num = DB.load_db_bases_to_num()
    all_bases_tuple, all_bases_num = db_bases_to_num.keys(), db_bases_to_num.values()
    del db_bases_to_num
    all_bases_tuple = [x for (y,x) in sorted(zip(all_bases_num,all_bases_tuple))]
    del all_bases_num

    for afile in onlyfiles:
        db = DB.load_db_pickle(join(path,afile))
        name = all_bases_tuple[db.values()[0].pop()]
        rp = rail_place.RailPlace(name)
        name = rp.getXYZ(turb)
        name = [round(name[0],3), round(name[1],3), round(name[2],3)]
        name = str(name)
        name = name.replace(', ','_')
        name = name.replace('[','')
        name = name.replace(']','')
        DB.save_db_pickle(db, join(path,name+'.pkl'))

def plot_robot_in_base():
    db = DB.load_db_bases_to_num()

def generate_robot_positions():
    rp = rail_place.rand_rail(turb, 1000)
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
    xml_trajectories_path = join(blade_folder,"trajectory/trajectory.xml")
    blade = blade_modeling.BladeModeling(turb, turb.blades[0])
    blade.load_trajectory(xml_trajectories_path)
    DB = db.DB(directory, blade)
    del blade
    return

if __name__ == '__main__':
    dir_test = os.path.join(os.path.realpath('.'),'test')
    os.environ['OPENRAVE_DATA'] = str(dir_test)
    cfg = TurbineConfig.load('turbine_unittest.cfg','test')
    turb = Turbine(cfg)
    #DB = db.DB('db')#, blade)
    vis = Visualizer(turb.env)

    #rename_files_by_base()
    #convert()
    #generate()
    #merge()



