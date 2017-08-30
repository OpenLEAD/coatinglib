from turbine import Turbine
from turbine_config import TurbineConfig
import os
from os.path import join, isfile, realpath
from numpy import array

dir_test = join(realpath('.'),'test')
os.environ['OPENRAVE_DATA'] = str(dir_test)
cfg = TurbineConfig.load('turbine_unittest.cfg','test')
turbine = Turbine(cfg)

import visualizer
vis = visualizer.Visualizer(turbine.env)

import blade_coverage
import db

grid = 77

robot = turbine.robot
manip = robot.GetActiveManipulator()
robot.GetLink('Flame').Enable(False)
DB = db.DB('FACE',turbine)

T = array([[ 1.,  0.,  0.,  0.],
       [ 0.,  0.,  1.,  0.],
       [ 0., -1.,  0.,  0.],
       [ 0.,  0.,  0.,  1.]])
DB.T = T

#psa = (1.8000000000000034, 0.10879140586688345, 1.349065850398866)

psa = (-2.0, -0.72551222845528596, 1.3490658503988662)


threshold = 5e-2
path = blade_coverage.base_grid_validation(turbine, psa, DB, grid, threshold)

