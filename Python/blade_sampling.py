from coating import *
from numpy import *
from openravepy import *
from openravepy.misc import SpaceSamplerExtra
import time
import math

env=Environment()
env.Load("/home/renan/workspace/coatinglib/Turbina/env_mh12_0_16.xml")
target = env.GetBodies()[0]

# PARAMETERS
delta = 0.001
normalanglerange = 0
directiondelta = 0.2
coatingdistance = 0.23

#MAIN
approachrays,_,_,_,_=PointsToReach(env, target, delta, normalanglerange,directiondelta,'s2')
approachrays=PointsForCoating(approachrays,coatingdistance)

savez_compressed('blade_sampling_full/s2.npz', array=approachrays)
