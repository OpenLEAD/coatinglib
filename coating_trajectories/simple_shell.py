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
import mathtools
grid = 0
robot = turb.robot
manip = robot.GetActiveManipulator()
robot.GetLink('Flame').Enable(False)
DB = db.DB('LIP',turb)

T = array([[ 1.,  0.,  0.,  0.],
       [ 0.,  0.,  1.,  0.],
       [ 0., -1.,  0.,  0.],
       [ 0.,  0.,  0.,  1.]])

psa = (1.0000000000000027, 1.0517976918168939, 0.30186829920226843)

for blade in turb.blades:
    blade.SetTransform(T)
    
rp = rail_place.RailPlace(psa)
turb.place_rail(rp)
turb.place_robot(rp)

threshold = 5e-2
organized_rays_list = blade_coverage.organize_rays_in_parallels(DB, grid)
joint_path, linear_interpolated_rays = blade_coverage.move_dijkstra(turb, DB.load_blade(), organized_rays_list, threshold)
joint_path_2 = blade_coverage.refine_dijkstra(turb, joint_path, linear_interpolated_rays, threshold)
new_joint_path = blade_coverage.smooth_joint_MLS(turb, joint_path_2)
for path in new_joint_path:
	for joint in path:
		robot.SetDOFValues(joint)
		time.sleep(0.05)
