from numpy import *
import implicit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

Rays1 = load('coated_points/pos4_black.npz')
Rays1  = Rays1['array']
Rays2 = load('coated_points/pos4_blue.npz')
Rays2  = Rays2['array']
rR = array(concatenate((Rays1,Rays2)))

def plot():
    ax = implicit.plot_implicit_load('implicit/QP/sigma_5/',bbox=(-4,4,-4,4,-4,4))
    Axes3D.scatter(ax,rR[:,0],rR[:,1],rR[:,2])
    plt.show()
