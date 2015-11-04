import coating
from openravepy import *
from numpy import *
from math import *
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import matplotlib.cm as cmx
from mpl_toolkits.mplot3d import Axes3D

# Get poinst for coating
velocities = load('VelList_fullHD1.npz')
velocities = velocities['array']

thetas = load('ThetaList_fullHD1.npz')
thetas = thetas['array']

omegas = load('OmegaList_fullHD1.npz')
omegas  = omegas['array']

nearlist = load('nearPointsByNumberOfPoints1_fullHD.npz')
nearlist  = nearlist['array']

manipulabilities = load('Manipulabilities1_fullHD.npz')
manipulabilities  = manipulabilities['array']

referencePoints = nearlist[:,0]
maxmanipulabilities = []
for i in range(0,len(manipulabilities)):
    maxmanipulabilities.append(max(manipulabilities[i]))

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

scatter3d(referencePoints[:,0],referencePoints[:,1],referencePoints[:,2], maxmanipulabilities, colorsMap='jet')
