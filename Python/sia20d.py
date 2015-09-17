import openravepy
from numpy import *
from openravepy import *
env=Environment()
env.SetViewer('qtcoin')
env.Load("/home/renan/Documents/EMMA/SIA20D/SIA20D.xml")
robot = env.GetRobots()[0]


with env:
    ikmodel = databases.inversekinematics.InverseKinematicsModel(robot=robot,iktype=IkParameterization.Type.Transform6D)
    if not ikmodel.load():
        ikmodel.autogenerate()
