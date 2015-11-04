import coating
from openravepy import *
from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from math import *

# Get poinst for coating
omegas = load('NewOmegas0_HD.npz')
omegas  = omegas['array']

alltriopoints = load('alltriopoints0_HD.npz')
alltriopoints  = alltriopoints['array']

x=[];y=[];z=[];nx=[];ny=[];nz=[]
maxOmegas = array([220,200,220,410,410,610])*pi/180

for i in range(0,len(omegas)):
    io = 1
    for j in range(0,len(omegas[i])):
        for k in range(0,len(omegas[i][j])):
            if (omegas[i][j][k]<maxOmegas).all()==True:
                x.append(alltriopoints[i][0][0][0])
                y.append(alltriopoints[i][0][0][1])
                z.append(alltriopoints[i][0][0][2])
                io = 0
                break
    if io:
        nx.append(alltriopoints[i][0][0][0])
        ny.append(alltriopoints[i][0][0][1])
        nz.append(alltriopoints[i][0][0][2])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

Axes3D.scatter(ax, nx, ny, nz,s=40,c='r')
Axes3D.scatter(ax, x, y, z,c='g')
ax.axis("off")

plt.show()
