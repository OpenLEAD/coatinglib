from turbine import Turbine
from turbine_config import TurbineConfig
import os
from os.path import join, isfile, realpath

dir_test = join(realpath('.'),'test')
os.environ['OPENRAVE_DATA'] = str(dir_test)
cfg = TurbineConfig.load('turbine_unittest.cfg','test')
turb = Turbine(cfg)

import visualizer
vis = visualizer.Visualizer(turb.env)

import blade_coverage
import db
import planning
from openravepy import interfaces, IkParameterization, RaveCreateTrajectory

turb.robot.GetLink('Flame').Enable(False)
DB = db.DB('FACE',turb)
joint_solutions = blade_coverage.trajectory_generation(DB, 7)


