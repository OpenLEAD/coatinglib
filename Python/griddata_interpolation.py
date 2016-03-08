from numpy import *
from scipy.interpolate import griddata
import matplotlib.pyplot as plt


Rays1 = load('blade_sampling/blade_crop_fast.npz')
Rays1  = Rays1['array']
Rays2 = load('blade_sampling/blade_crop_fast2.npz')
Rays2  = Rays2['array']
rR = concatenate((Rays1,Rays2))

#values = zeros((len(rR),1))
grid_x, grid_y = mgrid[min(rR[:,0]):max(rR[:,0]):5000j,
                               min(rR[:,1]):max(rR[:,1]):5000j]
#grid_x, grid_y = mgrid[0:1:100j, 0:1:200j]

S = griddata(rR[:,0:2],rR[:,2], (grid_x,grid_y),method='linear')

#plt.imshow(S.T, extent=(0,1,0,1), origin='lower')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
