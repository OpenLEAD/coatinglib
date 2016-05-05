from numpy import *
import implicit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

reachableRays = load('blade_sampling/blade_crop_fast.npz')
reachableRays  = reachableRays['array']
AllreachableRays = load('blade_sampling/blade_crop_fast2.npz')
AllreachableRays  = AllreachableRays['array']
rays = concatenate((reachableRays,AllreachableRays))

v = array([  9.45575516e+00,  -1.76521543e+00,  -6.61304926e+00,  -7.97589786e+00,
       -5.32902294e-01,   9.05160504e+00,  1.07514601e+01,   1.08974083e+01,
       7.05698039e-01,   2.35020903e-02,  -2.15901680e+00,   1.48090437e-01,
       3.62980652e-01,   7.69228097e-02,   1.43424673e-03,   2.04874422e+01,
      -1.12191974e+01,  -1.06890244e+01,  -6.35856103e+00,   1.52680741e+01,
       9.57623095e+00,   7.60893369e+00,  -5.48101561e-01,  -2.56629783e+00,
       1.50685017e-01,   2.06482087e+01,  -8.42742799e+00,  -2.33634067e+00,
       1.13519018e+01,   1.66986026e+00,   7.24618668e-01,   9.57373566e+00,
      -1.97827712e+00,   2.32742998e+00,   2.05620074e+00])

def vector4(x,y,z):
    return [1, z, z**2, z**3, z**4, y, y*z, y*z**2, y*z**3, y**2, y**2*z, y**2*z**2, y**3, y**3*z, y**4, x, x*z, x*z**2, x*z**3, x*y, x*y*z, x*y*z**2, x*y**2, x*y**2*z, x*y**3, x**2, x**2*z, x**2*z**2, x**2*y, x**2*y*z, x**2*y**2, x**3, x**3*z, x**3*y, x**4]

@vectorize
def fn4(x,y,z):
    return dot(v,vector4(x,y,z))

def plot():
    ax = implicit.plot_implicit(fn4,bbox=(-5,5,-5,5,-5,5),rc=150,ns=20,colors='r')
    Axes3D.scatter(ax,rays[:,0],rays[:,1],rays[:,2])
    plt.show()
