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

##reachableRays = load('blade_sampling_full/blade_crop_ufast.npz')
##reachableRays  = reachableRays['array']
reachableRays = load('Blade/jiraublade_points.npz')
reachableRays  = reachableRays['array']
rays = reachableRays

rayspos = copy.copy(rays)
eps = 0.001
for i in range(0,len(rayspos)):
    rayspos[i][0:3] = rayspos[i][0:3]+rayspos[i][3:6]*eps

raysExt = concatenate((rays,rayspos))
N=len(raysExt)
#savetxt("raysExt.csv", raysExt[:,0:3], delimiter=",")
    
w=[]
#w = load('w/x3_full_ufast.npz')
#w = w['array']
w = load('Blade/RBF/jiraublade_r3_w.npz')
w = w['array']

#savetxt("w.csv", w, delimiter=",")

handles2=[]
env2=None
index=-1
def set_handles(handles,env):
    global handles2, env2
    handles2 = handles
    env2 = env
    return

def makeTree(rR):
    Tree = KDTree(rR[:,0:3])
    return Tree

Tree = makeTree(rays)

def poly(ci):
    return poly1(ci)

def poly1(ci):
    return [1, ci[0], ci[1], ci[2]]

def phi(ci,cj):
    return x3(ci,cj)

def dphi(ci,cj):
    return dx3(ci,cj)

def x3(ci,cj):
    c = ci-cj
    c = sqrt(dot(c,c))
    return c*c*c

def dx3(ci,cj):
    c = ci-cj
    nc = 3*sqrt(dot(c,c))
    return array([nc*c[0],nc*c[1],nc*c[2]])

def phi_log(ci,cj):
    c = ci-cj
    c2 = dot(c,c)
    if c2==0:
        return 0
    c = sqrt(c2)
    return c2*log(c)

def dphi_log(ci,cj):
    c = ci-cj
    c2 = dot(c,c)
    if c2==0:
        return array([0,0,0])
    return array([c[0]*log(c2)+c[0],c[1]*log(c2)+c[1],c[2]*log(c2)+c[2]])

def gauss(ci,cj):
    e = 0.1
    c = ci-cj
    c2 = e**2*dot(c,c)
    return exp(-c2)

def dgauss(ci,cj):
    e = 0.1; e2 = e*e 
    k = -2*e2
    c = ci-cj
    c2 = e2*dot(c,c)
    oexp = exp(-c2)
    return array([k*c[0]*oexp,k*c[1]*oexp,k*c[2]*oexp])

def IQ(ci,cj):
    e = 0.001
    c = ci-cj
    c2 = e**2*dot(c,c)
    return 1/(1+c2)

def dIQ(ci,cj):
    e2 = 0.001**2; k = -2*e2
    c = ci-cj
    c2 = 1/((e2*dot(c,c)+1)**2)
    return array([k*c[0]*c2,k*c[1]*c2,k*c[2]*c2])

def f(c):
    f=0
    for i in range(0,N):
        f+=w[i]*phi(c,raysExt[i][0:3])
    return f

def f3(x,y,z):
    return f([x,y,z])

def df(c):
    df=array([0,0,0])
    for i in range(0,N):
        df=df+w[i]*dphi(c,raysExt[i][0:3])
    return df
        
def plot():
    ax = implicit.plot_implicit(f3,bbox=(-5,5,-5,5,-5,5))
    Axes3D.scatter(ax,rays[:,0],rays[:,1],rays[:,2])
    plt.show()
