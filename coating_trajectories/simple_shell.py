from turbine import Turbine
from turbine_config import TurbineConfig
import os
from os.path import join, isfile, realpath
from numpy import array
from time import time

dir_test = join(realpath('.'),'test')
os.environ['OPENRAVE_DATA'] = str(dir_test)
cfg = TurbineConfig.load('turbine_unittest.cfg','test')
turbine = Turbine(cfg)

import visualizer
vis = visualizer.Visualizer(turbine.env)

import blade_coverage
import db

grid = 7
robot = turbine.robot
manip = robot.GetActiveManipulator()
robot.GetLink('Flame').Enable(False)
DB = db.DB('FACE',turbine)

T = array([[ 1.,  0.,  0.,  0.],
           [ 0.,  0.25881905,  0.96592583,  0.],
           [ 0., -0.96592583,  0.25881905,  0.],
           [ 0.,  0.,  0.,  1.]])
DB.T = T

psa = (1.8500000000000032, 0.53725006964517363, 0.82546707480056725)

threshold = 5e-2
t = time()
path = blade_coverage.base_grid_validation(turbine, psa, DB, grid, threshold)
print time()-t
