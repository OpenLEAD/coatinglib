from numpy import *
from openravepy import *
from openravepy.misc import SpaceSamplerExtra
from coating import *
import time
env=Environment()
env.Load("/home/renan/Documents/EMMA/Turbina/env_motoman.xml")
robot = env.GetRobots()[0]
target = env.GetBodies()[0]
manip = robot.GetActiveManipulator()

# PARAMETERS
delta = 0.02
normalanglerange = 0
directiondelta = 0.4
anglerange = pi/6
coatingdistance = 0.23
robottobladedistance = 0.8
pN = numpy.array([0.005,-3.422,-0.035,0.693,-0.02603,-0.72047])
facevector = [0,0,-1]
initialdistance = 0.23
theta = [0,0,0]

#MAIN
with env:
    ikmodel = databases.inversekinematics.InverseKinematicsModel(robot=robot,
                                                                 iktype=IkParameterization.Type.Transform6D)
    if not ikmodel.load():
        ikmodel.autogenerate()

    approachrays = load('bladepoints.npz')
    approachrays = approachrays['array']
    N = approachrays.shape[0]
    Ttarget = target.GetTransform()

    reachableRays, bestDistance = BestBaseDistance(pN,initialdistance,approachrays,robot,ikmodel,facevector,theta) 
    print bestDistance
