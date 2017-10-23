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

#grid = 77
grid = 8
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

#psa = (1.8000000000000034, 1.5254013670473676, 1.0000000000000002)
psa = (1.7000000000000033, 0.73805000722783654, 0.82546707480056736)


#psa = (-2.0, -0.72551222845528596, 1.3490658503988662)


threshold = 5e-2
t = time()
path = blade_coverage.base_grid_validation(turbine, psa, DB, grid, threshold)
print time()-t
