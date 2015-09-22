import coating
from openravepy import *
from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from math import *

# Get poinst for coating
Torques = load('NewTorques0_HD.npz')
Torques = Torques['array']

alltriopoints = load('alltriopoints0_HD.npz')
alltriopoints  = alltriopoints['array']

maxTorques = array([9999,9999,9999,22,22,9.8])

nonFeasibleTPoints = []
nonFeasibleTTorques = []
FeasibleTPoints = []
FeasibleTTorques = []

for i in range(0,len(Torques)):
    io = 1
    for j in range(0,len(Torques[i])):
        if (Torques[i][j]<maxTorques).all()==True:
            FeasibleTPoints.append(alltriopoints[i])
            FeasibleTTorques.append(Torques[i])
            io = 0
            break
    if io:
        nonFeasibleTPoints.append(alltriopoints[i])
        nonFeasibleTTorques.append(Torques[i])
        
nonFeasibleTReferences = []
FeasibleTReferences = []
for nonFeasibleTPoint in nonFeasibleTPoints:
    nonFeasibleTReferences.append(nonFeasibleTPoint[0][0])
for FeasibleTPoint in FeasibleTPoints:
    FeasibleTReferences.append(FeasibleTPoint[0][0])    

nxarray = []
nyarray = []
nzarray = []
for nonFeasibleTReference in nonFeasibleTReferences:
    nxarray.append(nonFeasibleTReference[0])
    nyarray.append(nonFeasibleTReference[1])
    nzarray.append(nonFeasibleTReference[2])

xarray = []
yarray = []
zarray = []
for FeasibleTReference in FeasibleTReferences:
    xarray.append(FeasibleTReference[0])
    yarray.append(FeasibleTReference[1])
    zarray.append(FeasibleTReference[2])    

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

Axes3D.scatter(ax, nxarray, nyarray, nzarray,c='r')
#Axes3D.scatter(ax, xarray, yarray, zarray)
ax.axis("off")

plt.show()
