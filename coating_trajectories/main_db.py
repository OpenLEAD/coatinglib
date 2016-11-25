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
    folder = "jiraublade"
    xml_trajectories_path = os.path.join(folder,"trajectory/trajectory.xml")
    blade = blade_modeling.BladeModeling(turb, turb.blades[0])
    blade.load_trajectory(xml_trajectories_path)
    rpt = rail_place.RailPlace((.80,0.15,0))

    db_bases_to_num = DB.load_db_bases_to_num()
    
    while True:
        rp = rail_place.rand_rail(turb, 1)[0]
        turb.place_rail(rp)
        turb.place_robot(rp)

        if turb.check_rail_collision():
            continue
        if turb.check_robotbase_collision():
            continue

        break
        #for base in db_bases_to_num.keys():
            

    db = DB.generate_db(turb, blade, rp)
    now = datetime.now()
    DB.save_db_pickle(db,'db/servidor/nonconverted/' + now.strftime('%X').replace(':','_') + '.pkl')

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
        name = (round(name[0],3), round(name[1],3), round(name[2],3))
        name = str(name)
        name = name.replace(', ','_')
        name = name.replace('(','')
        name = name.replace(')','')
        DB.save_db_pickle(db, join(path,name+'.pkl'))

def plot_robot_in_base():
    db = DB.load_db_bases_to_num()

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



