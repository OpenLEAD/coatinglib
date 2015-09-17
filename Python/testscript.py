import coating
from openravepy import *
from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from math import *

# Get poinst for coating
nearlist = load('nearPointsByNumberOfPoints0_8.npz')
nearlist = nearlist['array']

allangles = load('allangles0_8.npz')
allangles = allangles['array']

omegas = load('omegas0.npz')
omegas = omegas['array']

NewOmegas = load('NewOmegas0.npz')
NewOmegas = NewOmegas['array']


alphas = load('alphas0.npz')
alphas = alphas['array']

alltriopoints = load('alltriopoints0.npz')
alltriopoints = alltriopoints['array']

thetas = load('thetas0.npz')
thetas = thetas['array']

maxOmega = 1.0*pi/180*array([220,200,220,410,410,610])

nonFeasiblePoints = []
nonFeasibleOmegas = []
nonFeasibleAlphas = []
nonFeasibleThetas = []
FeasiblePoints = []
FeasibleOmegas = []
FeasibleAlphas = []
FeasibleThetas = []


for i in range(0,len(omegas)):
    io = 1
    for j in range(0,len(omegas[i])):
        if (omegas[i][j][0]<maxOmega).all()==True and (omegas[i][j][1]<maxOmega).all()==True:
            FeasiblePoints.append(alltriopoints[i])
            FeasibleOmegas.append(omegas[i][j])
            FeasibleAlphas.append(alphas[i][j])
            FeasibleThetas.append(thetas[i][j])
            io = 0
            break
    if io:
        nonFeasiblePoints.append(alltriopoints[i])
        nonFeasibleOmegas.append(omegas[i])
        nonFeasibleAlphas.append(alphas[i])
        nonFeasibleThetas.append(thetas[i])

nonFeasibleReferences = []
FeasibleReferences = []
for nonFeasiblePoint in nonFeasiblePoints:
    nonFeasibleReferences.append(nonFeasiblePoint[0][0])
for FeasiblePoint in FeasiblePoints:
    FeasibleReferences.append(FeasiblePoint[0][0])    

nxarray = []
nyarray = []
nzarray = []
for nonFeasibleReference in nonFeasibleReferences:
    nxarray.append(nonFeasibleReference[0])
    nyarray.append(nonFeasibleReference[1])
    nzarray.append(nonFeasibleReference[2])

xarray = []
yarray = []
zarray = []
for FeasibleReference in FeasibleReferences:
    xarray.append(FeasibleReference[0])
    yarray.append(FeasibleReference[1])
    zarray.append(FeasibleReference[2])    

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

Axes3D.scatter(ax, nxarray, nyarray, nzarray,c='r')
#Axes3D.scatter(ax, xarray, yarray, zarray)
ax.axis("off")

#plt.show()
