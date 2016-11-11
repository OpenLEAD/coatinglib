#!/usr/bin/env python
import db
import blade_modeling
from turbine import Turbine
from turbine_config import TurbineConfig, ConfigFileError
import os
from visualizer import Visualizer
from os.path import join, isfile
from os import listdir
import rail_place
from numpy import array, random
from datetime import datetime


def convert():
    DB.convert_db_point_base_directory('db/servidor/nonconverted', 'db/converted')

def merge():
    DB.merge_db_directory('db/converted')

def plot_gradient():
    vis = Visualizer(turb.env)
    DB.plot_points_gradient(vis)
    
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

if __name__ == '__main__':
    dir_test = os.path.join(os.path.realpath('.'),'test')
    os.environ['OPENRAVE_DATA'] = str(dir_test)
    cfg = TurbineConfig.load('turbine_unittest.cfg','test')
    turb = Turbine(cfg)
    DB = db.DB('db')#, blade)
    generate()



