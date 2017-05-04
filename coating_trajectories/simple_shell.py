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
manip = turb.robot.GetActiveManipulator()
#basemanip = interfaces.BaseManipulation(turb.robot)

turb.robot.GetLink('Flame').Enable(False)
DB = db.DB('FACE',turb)
score, joint_solutions = blade_coverage.base_grid_validation(DB, 7)
traj = blade_coverage.movetohandposition(DB.turb.robot,joint_solutions)
#joints = []
#for i in range(0,traj.GetNumWaypoints()):
#    joints.append(traj.GetWaypoint(i)[0:6])
#points = []
#for joint in joints:
#    turb.robot.SetDOFValues(joint)
#    points.append(manip.GetTransform()[0:3,3])
#p = vis.plot(points,'p',(1,0,0))

for listT in traj:
    for T in listT:
        p = vis.plot(T[0:3,3],'p',(1,0,0))
        p = vis.plot_normal(list(T[0:3,3])+list(T[0:3,0]),'p',(0,0,1))


