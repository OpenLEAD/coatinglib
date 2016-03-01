import coating
from numpy import *
from openravepy import *
from openravepy.misc import SpaceSamplerExtra
import time

env=Environment()
env.Load("/home/renan/workspace/coatinglib/Turbina/env_mh12_0_16.xml")
target = env.GetBodies()[0]

# PARAMETERS
delta = 0.02
normalanglerange = 0
directiondelta = 0.4
anglerange = pi/6


#MAIN
with env:

    approachrays = load('blade_sampling/blade2.npz')
    approachrays = approachrays['array']
    N = approachrays.shape[0]
    Ttarget = target.GetTransform()

    R0 = 1.43**2
    Rf = 3.85**2

    approachrays2 = []
    for ray in approachrays:
        r = dot(ray[0:3],transpose(ray[0:3]))
        if r>R0 and r<Rf:
            approachrays2.append(ray)
    
    savez_compressed('blade_sampling/blade_crop2.npz', array=approachrays2)
