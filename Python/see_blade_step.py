import coating
from numpy import *
from openravepy import *
from openravepy.misc import SpaceSamplerExtra
import time
import math

env=Environment()
env.SetViewer('qtcoin')
env.Load("/home/renan/workspace/coatinglib/Turbina/env_mh12_0_16.xml")
target = env.GetBodies()[0]

# PARAMETERS
delta = 0.02
normalanglerange = 0
directiondelta = 0.4
anglerange = pi/6


#MAIN
with env:

    approachrays = load('blade_sampling/blade_pcl_crop.npz')
    approachrays = approachrays['array']
    N = approachrays.shape[0]
    Ttarget = target.GetTransform()

   
    gapproachrays = c_[dot(approachrays[0:N,0:3],transpose(Ttarget[0:3,0:3]))+tile(Ttarget[0:3,3],(N,1)),dot(approachrays[0:N,3:6],transpose(Ttarget[0:3,0:3]))]
handles=[]
for ray in gapproachrays:
    handles.append(env.plot3(points=ray[0:3],pointsize=5,colors=array((1,0,0))))
    time.sleep(0.05)
