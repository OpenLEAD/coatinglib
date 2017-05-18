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
import planning
from openravepy import interfaces, IkParameterization, RaveCreateTrajectory
import time
import rail_place
grid = 0
robot = turb.robot
manip = robot.GetActiveManipulator()
robot.GetLink('Flame').Enable(False)
DB = db.DB('LIP',turb)

T = array([[ 1.,  0.,  0.,  0.],
       [ 0.,  0.,  1.,  0.],
       [ 0., -1.,  0.,  0.],
       [ 0.,  0.,  0.,  1.]])

psa = (1.2000000000000028, 1.2154441840549897, 0.4764012244017013)

for blade in turb.blades:
    blade.SetTransform(T)
    
rp = rail_place.RailPlace(psa)
turb.place_rail(rp)
turb.place_robot(rp)

organized_rays_list = blade_coverage.organize_rays_in_parallels(DB, grid)
joint_path = blade_coverage.move_dijkstra(turb, DB.load_blade(), organized_rays_list, 3e-2)
for path in joint_path:
	for joint in path:
		robot.SetDOFValues(joint)
		time.sleep(0.05)
