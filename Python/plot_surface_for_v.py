from numpy import *
import implicit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

reachableRays = load('blade_sampling/blade_crop_fast.npz')
reachableRays  = reachableRays['array']
AllreachableRays = load('blade_sampling/blade_crop_fast2.npz')
AllreachableRays  = AllreachableRays['array']
rays = concatenate((reachableRays,AllreachableRays))

v = array([ -3.83453602e+04,   1.00354987e+05, -9.75943644e+04,   4.18957355e+04,
  -6.70585470e+03,  -4.05112645e+04,  7.95287655e+04,  -5.15816715e+04,
   1.10823557e+04,  -1.60205691e+04 ,  2.09743085e+04 , -6.80531378e+03,
  -2.80976042e+03,   1.84028353e+03,  -1.84372023e+02 , -1.43173084e+05,
   2.79789333e+05,  -1.80792493e+05 ,  3.86818959e+04 , -1.13652775e+05,
   1.48085692e+05 , -4.78744151e+04 , -3.00276035e+04 ,  1.95681517e+04,
  -2.63928135e+03 , -1.99747945e+05 ,  2.59239434e+05 , -8.34693306e+04,
  -1.05878597e+05,   6.87244827e+04 , -1.40150336e+04 , -1.23469085e+05,
   7.98234087e+04 , -3.27683167e+04,  -2.85359378e+04])

def vector4(x,y,z):
    return [1, z, z**2, z**3, z**4, y, y*z, y*z**2, y*z**3, y**2, y**2*z, y**2*z**2, y**3, y**3*z, y**4, x, x*z, x*z**2, x*z**3, x*y, x*y*z, x*y*z**2, x*y**2, x*y**2*z, x*y**3, x**2, x**2*z, x**2*z**2, x**2*y, x**2*y*z, x**2*y**2, x**3, x**3*z, x**3*y, x**4]

@vectorize
def fn4(x,y,z):
    return dot(v,vector4(x,y,z))

def plot():
    ax = implicit.plot_implicit(fn4,bbox=(-5,5,-5,5,-5,5),rc=150,ns=50,colors='r')
    Axes3D.scatter(ax,rays[:,0],rays[:,1],rays[:,2])
    plt.show()
