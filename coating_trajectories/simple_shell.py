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

psa = (1.9000000000000035, 1.1917548675411605, 0.95993108859688125)

for blade in turb.blades:
    blade.SetTransform(T)
    
rp = rail_place.RailPlace(psa)
turb.place_rail(rp)
turb.place_robot(rp)



##score, joint_solutions_list = blade_coverage.base_grid_validation_parallel(DB, 7)
##traj = blade_coverage.movetohandposition_parallels(DB.turb.robot,joint_solutions)
##joints = []
##for i in range(0,traj.GetNumWaypoints()):
##    joints.append(traj.GetWaypoint(i)[7:7+6])
##points = []
##for joint in joints:
##    turb.robot.SetDOFValues(joint)
##    points.append(manip.GetTransform()[0:3,3])
##p = vis.plot(points,'p',(1,0,0))

organized_rays = blade_coverage.organize_rays_in_parallels(DB, grid)
joint_path = blade_coverage.move_dijkstra(turb, DB.load_blade(), organized_rays)
for path in joint_path:
	for joint in path:
		robot.SetDOFValues(joint)
		time.sleep(0.01)
