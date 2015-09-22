import coating
from openravepy import *
from numpy import *
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import matplotlib.cm as cmx
from mpl_toolkits.mplot3d import Axes3D
from math import *

Torques = load('Torques0_HD.npz')
Torques = Torques['array']

alltriopoints = load('alltriopoints0_HD.npz')
alltriopoints  = alltriopoints['array']

Max = []
nMax = []
xarray = []
yarray = []
zarray = []
xnarray = []
ynarray = []
znarray = []

for i in range(0,len(Torques)):
    t=[]
    for j in range(0,len(Torques[i])):
        t.append(array([Torques[i][j][3],Torques[i][j][4]]))
    io=1    
    for j in range(0,len(t)):
        if t[j][0]<=22 and t[j][1]<=22:
            Max.append(amax(t[j]))    
            xarray.append(alltriopoints[i][0][0][0])
            yarray.append(alltriopoints[i][0][0][1])
            zarray.append(alltriopoints[i][0][0][2])
            io=0
            break
    if io:    
        nMax.append(amin(t))      
        xnarray.append(alltriopoints[i][0][0][0])
        ynarray.append(alltriopoints[i][0][0][1])
        znarray.append(alltriopoints[i][0][0][2])
        
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

Axes3D.scatter(ax, xnarray, ynarray, znarray,c='r')
Axes3D.scatter(ax, xarray, yarray, zarray)
ax.axis("off")
plt.show()

def scatter3d(x,y,z, cs, colorsMap='jet'):
    cm = plt.get_cmap(colorsMap)
    cNorm = clr.Normalize(vmin=min(cs), vmax=max(cs))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(x, y, z, c=scalarMap.to_rgba(cs))
    scalarMap.set_array(cs)
    fig.colorbar(scalarMap)
    ax.axis("off")
    plt.show()

#scatter3d(xnarray,ynarray,znarray, nMax, colorsMap='jet')
