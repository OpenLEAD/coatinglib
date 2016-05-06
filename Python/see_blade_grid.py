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

    approachrays = load('blade_sampling_full/blade_crop_fast.npz')
    approachrays = approachrays['array']
    N = approachrays.shape[0]
    Ttarget = target.GetTransform()

    #Blade PLC:5
##    p = [0,0,-1.265741]
##    Rt = dot(dot(coating.RunitY(-0.5*math.pi),coating.RunitZ(0.5*math.pi)),coating.RunitY(-31*math.pi/180))
##    Ttarget = array([[Rt[0][0],Rt[0][1],Rt[0][2],p[0]],
##                     [Rt[1][0],Rt[1][1],Rt[1][2],p[1]],
##                     [Rt[2][0],Rt[2][1],Rt[2][2],p[2]],
##                     [0,0,0,1]])
##    
##    Rx=coating.RunitX(-0.5*math.pi)
##    T = array([[Rx[0][0],Rx[0][1],Rx[0][2],0],
##                     [Rx[1][0],Rx[1][1],Rx[1][2],0],
##                     [Rx[2][0],Rx[2][1],Rx[2][2],0],
##                     [0,0,0,1]])
##
##    Ttarget=dot(T,Ttarget)
    #---------------
    
    gapproachrays = c_[dot(approachrays[0:N,0:3],transpose(Ttarget[0:3,0:3]))+tile(Ttarget[0:3,3],(N,1)),dot(approachrays[0:N,3:6],transpose(Ttarget[0:3,0:3]))]
##    savez_compressed('blade_sampling/blade_pcl.npz', array=gapproachrays)
    approachgraphs = env.plot3(points=gapproachrays[:,0:3],pointsize=5,colors=array((1,0,0)))
