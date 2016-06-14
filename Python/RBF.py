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

reachableRays = load('blade_sampling/blade_crop_fast.npz')
#AllreachableRays = load('blade_sampling/blade_crop_fast2.npz')
reachableRays  = reachableRays['array']
#AllreachableRays  = AllreachableRays['array']
#rays = concatenate((reachableRays,AllreachableRays))
rays = reachableRays

#rays = load('blade_sampling_full/blade_crop_fast.npz')
#rays  = rays['array']
rayspos = copy.copy(rays)
raysneg = copy.copy(rays)
eps = 0.001
for i in range(0,len(rayspos)):
    rayspos[i][0:3] = rayspos[i][0:3]+rayspos[i][3:6]*eps
    raysneg[i][0:3] = raysneg[i][0:3]-raysneg[i][3:6]*eps
raysE = concatenate((rayspos,raysneg))

raysExt = concatenate((rays,raysE))
N=len(raysExt)
#savetxt("raysExt.csv", raysExt[:,0:3], delimiter=",")
    
w=[]
w = load('w/r3eps.npz')
w = w['array']
savetxt("w.csv", w, delimiter=",")

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

def update_w():
    global w
    w = compute_w()
    return
    
def poly(ci):
    return poly1(ci)

def poly1(ci):
    return [1, ci[0], ci[1], ci[2]]

def phi(ci,cj):
    return phi_x3(ci,cj)

def dphi(ci,cj):
    return dphi_x3(ci,cj)

def phi_x3(ci,cj):
    c = ci-cj
    c = sqrt(dot(c,c))
    return c**3

def dphi_x3(ci,cj):
    c = ci-cj
    nc = 3*sqrt(dot(c,c))
    return array([nc*c[0],nc*c[1],nc*c[2]])

def compute_w():
    K = []
    d = []
    for i in range(0,N):
        ki=[]
        for ray in raysExt:
            ki.append(phi(rays[i][0:3],ray[0:3]))
        K.append(ki)    
        if i>=N/2:
            d.append(eps)
        else: d.append(0)    
    return solve(K,d)    

def main():
    update_w()

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

def optmizeTan(p0, pnew, tol):
    def func(P):
        a=P-pnew
        return dot(a,a)
    def func_deriv(P):
        return 2*(P-pnew)
    def consfunc(P):
        return f(P)
    def consfunc_deriv(P):
        return df(P)
    cons = ({'type':'eq',
             'fun': consfunc,
            'jac':consfunc_deriv})
    res = minimize(func, p0,jac=func_deriv,constraints=cons, method='SLSQP',tol=tol, options={'disp': False})
    return res      
