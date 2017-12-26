import openravepy
from numpy import *
from openravepy import *

env = Environment()
env.Load('arm_example.xml')

robot = env.GetRobots()[0]
robot.Enable(False)

reachmodel = databases.kinematicreachability.ReachabilityModel(robot)

