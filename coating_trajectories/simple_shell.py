from turbine import Turbine
from turbine_config import TurbineConfig
import os
from os.path import join, isfile, realpath
from numpy import array

dir_test = join(realpath('.'),'test')
os.environ['OPENRAVE_DATA'] = str(dir_test)
cfg = TurbineConfig.load('turbine_unittest.cfg','test')
turb = Turbine(cfg)

import visualizer
vis = visualizer.Visualizer(turb.env)

import blade_coverage
import db

grid = 1
robot = turb.robot
manip = robot.GetActiveManipulator()
robot.GetLink('Flame').Enable(False)
DB = db.DB('LIP',turb)

T = array([[ 1.,  0.,  0.,  0.],
       [ 0.,  0.,  1.,  0.],
       [ 0., -1.,  0.,  0.],
       [ 0.,  0.,  0.,  1.]])
DB.T = T
psa = (1.7500000000000036, 0.55169264516096073, 1.3490658503988662)

threshold = 5e-2
path = blade_coverage.base_grid_validation(turb, psa, DB, grid, threshold)

