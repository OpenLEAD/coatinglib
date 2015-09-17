from numpy import *
from openravepy import *
from openravepy.misc import SpaceSamplerExtra
from coating import *
import time
env=Environment()
env.Load("/home/renan/Documents/EMMA/Turbina/env_SIA20D.xml")
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
#pN = numpy.array([-0.00912,-2.72560,-0.07085,0.66103,-0.02598,-0.74991])
#pN = numpy.array([0.005,-3.422,-0.035,0.693,-0.02603,-0.72047])
pN = numpy.array([-0.03656,-2.98442,-0.08606,0.66103,-0.02598,-0.74991])
initialdistance = 0.23
facevector = [0,1,0]
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
