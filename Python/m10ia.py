import openravepy
from numpy import *
from openravepy import *
env=Environment()
env.SetViewer('qtcoin')
env.Load("/home/renan/Documents/EMMA/M-10iA/M10iA.xml")
robot = env.GetRobots()[0]
manip = robot.GetActiveManipulator()

with env:
    ikmodel = databases.inversekinematics.InverseKinematicsModel(robot=robot,
                                                                 iktype=IkParameterization.Type.Transform6D)#, freeindices=[0])
    if not ikmodel.load():
        ikmodel.autogenerate()
