import coating
from openravepy import *
from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from math import *

# Get poinst for coating
Torques = load('Torques0_HD.npz')
Torques = Torques['array']

FeasiblePoints = load('FeasiblePoints0_HD.npz')
FeasiblePoints  = FeasiblePoints['array']

FeasibleOmegas = load('FeasibleOmegas0_HD.npz')
FeasibleOmegas  = FeasibleOmegas['array']

FeasibleAlphas = load('FeasibleAlphas0_HD.npz')
FeasibleAlphas  = FeasibleAlphas['array']

FeasibleThetas = load('FeasibleThetas0_HD.npz')
FeasibleThetas = FeasibleThetas['array']

maxTorques = array([9999,9999,9999,22,22,9.8])

nonFeasibleTPoints = []
nonFeasibleTOmegas = []
nonFeasibleTAlphas = []
nonFeasibleTThetas = []
FeasibleTPoints = []
FeasibleTOmegas = []
FeasibleTAlphas = []
FeasibleTThetas = []


for i in range(0,len(Torques)):
    io = 1
    for j in range(0,len(Torques[i])):
        if (Torques[i][j]<maxTorques).all()==True:
            FeasibleTPoints.append(FeasiblePoints[i])
            FeasibleTOmegas.append(FeasibleOmegas[i])
            FeasibleTAlphas.append(FeasibleAlphas[i])
            FeasibleTThetas.append(FeasibleThetas[i])
            io = 0
            break
    if io:
        nonFeasibleTPoints.append(FeasiblePoints[i])
        nonFeasibleTOmegas.append(FeasibleOmegas[i])
        nonFeasibleTAlphas.append(FeasibleAlphas[i])
        nonFeasibleTThetas.append(FeasibleThetas[i])

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
Axes3D.scatter(ax, xarray, yarray, zarray)
ax.axis("off")

savez_compressed('FeasibleTPoints0_HD.npz', array=FeasibleTPoints)
savez_compressed('FeasibleTAlphas0_HD.npz', array=FeasibleTAlphas)
savez_compressed('FeasibleTOmegas0_HD.npz', array=FeasibleTOmegas)
savez_compressed('FeasibleTThetas0_HD.npz', array=FeasibleTThetas)

plt.show()
