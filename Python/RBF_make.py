from numpy import *
from math import pi
import copy
from scipy.linalg import solve
import implicit
import time
from scipy.spatial import KDTree
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import minimize

reachableRays = load('blade_sampling_full/blade_crop_ufast.npz')
reachableRays  = reachableRays['array']
rays = reachableRays
rayspos = copy.copy(rays)
eps = 0.001
for i in range(0,len(rayspos)):
    rayspos[i][0:3] = rayspos[i][0:3]+rayspos[i][3:6]*eps
rays = concatenate((rays,rayspos))
N=len(rays)
 
def update_w():
    return compute_w()
    
def poly(ci):
    return poly1(ci)

def poly1(ci):
    return [1, ci[0], ci[1], ci[2]]

def phi(ci,cj):
    return x3(ci,cj)

def x3(ci,cj):
    c = ci-cj
    c = sqrt(dot(c,c))
    return c**3

def phi_log(ci,cj):
    c = ci-cj
    c2 = dot(c,c)
    if c2==0:
        return 0
    c = sqrt(c2)
    return c2*log(c)

def gauss(ci,cj):
    e = 0.1
    c = ci-cj
    c2 = e**2*dot(c,c)
    return exp(-c2)

def IQ(ci,cj):
    e = 0.001
    c = ci-cj
    c2 = e**2*dot(c,c)
    return 1/(1+c2)

def compute_w():
    K = []
    d = []
    for i in range(0,N):
        ki=[]
        for ray in rays:
            ki.append(phi(rays[i][0:3],ray[0:3]))
        K.append(ki)    
        if i>=N/2:
            d.append(eps)
        else: d.append(0)    
    return solve(K,d)    

def main():
    w = update_w()
    savez_compressed('w/'+'x3_full_ufast.npz', array=w)
    return

def f(c):
    f=0
    for i in range(0,N):
        f+=w[i]*phi(c,rays[i][0:3])
    return f

def f3(x,y,z):
    return f([x,y,z])
        
def plot():
    ax = implicit.plot_implicit(f3,bbox=(-5,5,-5,5,-5,5))
    Axes3D.scatter(ax,rays[:,0],rays[:,1],rays[:,2])
    plt.show()

if __name__ == '__main__':
    main()
