from turbine import Turbine
from turbine_config import TurbineConfig
import os
from os.path import join, isfile, realpath
from numpy import array
from time import time

cfg = TurbineConfig.load('test/turbine_unittest.cfg')
turbine = Turbine(cfg)

import visualizer
vis = visualizer.Visualizer(turbine.env)

import blade_coverage
import db

grid = 7
robot = turbine.robot
manip = robot.GetActiveManipulator()
robot.GetLink('Flame').Enable(False)
DB = db.DB(join(os.environ['PYTHON_DATABASE'],'FACE'),turbine)

T = turbine.blades[3].GetTransform()

DB.T = T
psa = (1.8500000000000032, 0.53725006964517363, 0.82546707480056725)

threshold = 5e-2
t = time()
path = blade_coverage.base_grid_validation(turbine, psa, DB, grid, threshold)
print time()-t
