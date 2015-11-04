import coating
from numpy import *
from openravepy import *
from openravepy.misc import SpaceSamplerExtra
import time

approachrays = load('blade2.npz')
approachrays = approachrays['array']

i=j=0
while True:
    while True:
        dist = linalg.norm(array([approachrays[i][0],approachrays[i][1],approachrays[i][2]])-array([approachrays[j][0],approachrays[j][1],approachrays[j][2]]))
        if dist<0.01:
            approachrays = delete(approachrays,j,0)
        j+=1
        if j==len(approachrays):break
    i+=1
    j=0
    if i==len(approachrays):break
    print str(i)+'/'+str(len(approachrays))
            
