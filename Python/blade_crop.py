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

    r = load('blade_sampling_full/r.npz')
    r = r['array']
    approachrays = r


    N = approachrays.shape[0]
    Ttarget = target.GetTransform()

    approachrays3=[]
    for ray in approachrays:
        approachrays3.append([ray[0]+ray[3]*-0.23,ray[1]+ray[4]*-0.23,ray[2]+ray[5]*-0.23,ray[3],ray[4],ray[5]])

    approachrays3=array(approachrays3)
    R0 = 1.59**2
    Rf = 3.75**2

    approachrays2 = []
    for ray in approachrays3:
        r = dot(ray[0:3],transpose(ray[0:3]))
        if r>R0 and r<Rf:
            approachrays2.append(ray)

    approachrays3=[]
    for ray in approachrays2:
        approachrays3.append([ray[0]+ray[3]*+0.23,ray[1]+ray[4]*+0.23,ray[2]+ray[5]*+0.23,ray[3],ray[4],ray[5]])
    
    savez_compressed('blade_sampling_full/r_crop.npz', array=approachrays3)
