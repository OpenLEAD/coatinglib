import openravepy
from numpy import *
from openravepy import *
env=Environment()
env.SetViewer('qtcoin')
env.Load("/home/renan/Documents/EMMA/motoman/motoman_mh122.xml")
robot = env.GetRobots()[0]
RaveSetDebugLevel(DebugLevel.Debug)

with env:
    ikmodel = databases.inversekinematics.InverseKinematicsModel(robot=robot,iktype=IkParameterization.Type.Transform6D)
    if not ikmodel.load():
        ikmodel.autogenerate()

