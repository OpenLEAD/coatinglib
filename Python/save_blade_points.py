from coating import *
from numpy import *
from openravepy import *
from openravepy.misc import SpaceSamplerExtra
import time

env=Environment()
env.Load("/home/renan/Documents/EMMA/Turbina/env_motoman0.xml")
target = env.GetBodies()[0]

# PARAMETERS
delta = 0.02
normalanglerange = 0
directiondelta = 0.4
coatingdistance = 0.23

#MAIN
approachrays=PointsToReach(env, target, delta, normalanglerange,directiondelta,'r')
approachrays=PointsForCoating(approachrays,coatingdistance)

savez_compressed('bladepointsMotoman0.npz', array=approachrays)
