from coating import *
from numpy import *
from openravepy import *
from openravepy.misc import SpaceSamplerExtra
import time

env=Environment()
env.Load("/home/renan/Documents/EMMA/Turbina/env_mh12_0_16.xml")
target = env.GetBodies()[0]

# PARAMETERS
delta = 0.01
normalanglerange = 0
directiondelta = 0.4
coatingdistance = 0.23

#MAIN
approachrays=PointsToReach(env, target, delta, normalanglerange,directiondelta,'b')
approachrays=PointsForCoating(approachrays,coatingdistance)

savez_compressed('bladepointsTEST2.npz', array=approachrays)
