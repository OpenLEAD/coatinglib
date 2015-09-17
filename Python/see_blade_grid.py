import coating
from numpy import *
from openravepy import *
from openravepy.misc import SpaceSamplerExtra
import time

env=Environment()
env.SetViewer('qtcoin')
env.Load("/home/renan/Documents/EMMA/Turbina/env_mh12_0_16.xml")
#env.Load("/home/renan/Documents/EMMA/Turbina/env_motoman0.xml")
target = env.GetBodies()[0]

# PARAMETERS
delta = 0.02
normalanglerange = 0
directiondelta = 0.4
anglerange = pi/6


#MAIN
with env:
   # approachrays = load('bladepointslr2.npz')
   # approachrays = approachrays['array']
   # N = approachrays.shape[0]
   # Ttarget = target.GetTransform()

    #approachrays2 = load('bladepointsLeft2.npz')
    approachrays2 = load('bladepointsTEST2.npz')
    approachrays2 = approachrays2['array']
    N2 = approachrays2.shape[0]
    Ttarget = target.GetTransform()
    
    #PLOT BLADE POINTS FOR COATING
    #gapproachrays = c_[dot(approachrays[0:N,0:3],transpose(Ttarget[0:3,0:3]))+tile(Ttarget[0:3,3],(N,1)),dot(approachrays[0:N,3:6],transpose(Ttarget[0:3,0:3]))]
    #approachgraphs = env.plot3(points=gapproachrays[:,0:3],pointsize=5,colors=array((1,0,0)))

    gapproachrays2 = c_[dot(approachrays2[0:N2,0:3],transpose(Ttarget[0:3,0:3]))+tile(Ttarget[0:3,3],(N2,1)),dot(approachrays2[0:N2,3:6],transpose(Ttarget[0:3,0:3]))]
    approachgraphs2 = env.plot3(points=gapproachrays2[:,0:3],pointsize=5,colors=array((1,0,0)))
