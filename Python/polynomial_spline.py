from numpy import *
from math import pi
import copy
from scipy.linalg import solve
import implicit
import time
from scipy.spatial import KDTree
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

##reachableRays = load('coated_points/pos4_black.npz')
##reachableRays  = reachableRays['array']
##AllreachableRays = load('coated_points/pos4_blue.npz')
##AllreachableRays  = AllreachableRays['array']
##rays = concatenate((reachableRays,AllreachableRays))

reachableRays = load('blade_sampling/blade_crop_fast.npz')
reachableRays  = reachableRays['array']
AllreachableRays = load('blade_sampling/blade_crop_fast2.npz')
AllreachableRays  = AllreachableRays['array']
rays = concatenate((reachableRays,AllreachableRays))


v=[]
r=0.4
center=[500,0,0]
idx=[]


def load_rays(lrays):
    global rays
    rays=lrays
    return

def update_v(newv):
    global v
    v=newv
    return

def update_idx(newidx):
    global idx
    idx=newidx
    return

def update_center(newcenter):
    global center
    center=newcenter
    return

def makeTree(rR):
    Tree = KDTree(rR[:,0:3])
    return Tree

Tree = makeTree(rays)

def getSamples(point,r):
    idx = Tree.query_ball_point(point,r)
    return idx

def vector4(x,y,z):
    return [1, z, z**2, z**3, z**4, y, y*z, y*z**2, y*z**3, y**2, y**2*z, y**2*z**2, y**3, y**3*z, y**4, x, x*z, x*z**2, x*z**3, x*y, x*y*z, x*y*z**2, x*y**2, x*y**2*z, x*y**3, x**2, x**2*z, x**2*z**2, x**2*y, x**2*y*z, x**2*y**2, x**3, x**3*z, x**3*y, x**4]

def dvector4(x,y,z):
    dx = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, z, z**2, z**3, y, y*z, y*z**2, y**2, y**2*z, y**3, 2*x, 2*x*z, 2*x*z**2, 2*x*y, 2*x*y*z, 2*x*y**2, 3*x**2, 3*x**2*z, 3*x**2*y, 4*x**3]
    dy = [0, 0, 0, 0, 0, 1, z, z**2, z**3, 2*y, 2*y*z, 2*y*z**2, 3*y**2, 3*y**2*z, 4*y**3, 0, 0, 0, 0, x, x*z, x*z**2, 2*x*y, 2*x*y*z, 3*x*y**2, 0, 0, 0, x**2, x**2*z, 2*x**2*y, 0, 0, x**3, 0]
    dz = [0, 1, 2*z, 3*z**2, 4*z**3, 0, y, 2*y*z, 3*y*z**2, 0, y**2, 2*y**2*z, 0, y**3, 0, 0, x, 2*x*z, 3*x*z**2, 0, x*y, 2*x*y*z, 0, x*y**2, 0, 0, x**2, 2*x**2*z, 0, x**2*y, 0, 0, x**3, 0, 0]
    G=[dx,dy,dz]
    return G

def ddvector4(x,y,z):
    dxdx = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2*z, 2*z**2, 2*y, 2*y*z, 2*y**2, 6*x, 6*x*z, 6*x*y, 12*x**2]
    dxdy = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, z, z**2, 2*y, 2*y*z, 3*y**2, 0, 0, 0, 2*x, 2*x*z, 4*x*y, 0, 0, 3*x**2, 0]
    dxdz = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2*z, 3*z**2, 0, y, 2*y*z, 0, y**2, 0, 0, 2*x, 4*x*z, 0, 2*x*y, 0, 0, 3*x**2, 0, 0]
    #dydx = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, z, z**2, 2*y, 2*y*z, 3*y**2, 0, 0, 0, 2*x, 2*x*z, 4*x*y, 0, 0, 3*x**2, 0]
    dydy = [0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2*z, 2*z**2, 6*y, 6*y*z, 12*y**2, 0, 0, 0, 0, 0, 0, 0, 2*x, 2*x*z, 6*x*y, 0, 0, 0, 0, 0, 2*x**2, 0, 0, 0, 0]
    dydz = [0, 0, 0, 0, 0, 0, 1, 2*z, 3*z**2, 0, 2*y, 4*y*z, 0, 3*y**2, 0, 0, 0, 0, 0, 0, x, 2*x*z, 0, 2*x*y, 0, 0, 0, 0, 0, x**2, 0, 0, 0, 0, 0]
    #dzdx = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2*z, 3*z**2, 0, y, 2*y*z, 0, y**2, 0, 0, 2*x, 4*x*z, 0, 2*x*y, 0, 0, 3*x**2, 0, 0]
    #dzdy = [0, 0, 0, 0, 0, 0, 1, 2*z, 3*z**2, 0, 2*y, 4*y*z, 0, 3*y**2, 0, 0, 0, 0, 0, 0, x, 2*x*z, 0, 2*x*y, 0, 0, 0, 0, 0, x**2, 0, 0, 0, 0, 0]
    dzdz = [0, 0, 2, 6*z, 12*z**2, 0, 0, 2*y, 6*y*z, 0, 0, 2*y**2, 0, 0, 0, 0, 0, 2*x, 6*x*z, 0, 0, 2*x*y, 0, 0, 0, 0, 0, 2*x**2, 0, 0, 0, 0, 0, 0, 0]
    return [dxdx,dxdy,dxdz,dydy,dydz,dzdz]

def matrix():
    m = []
    S = []
    for ind in idx:
        a=vector4(rays[ind][0],rays[ind][1],rays[ind][2])   
        da=dvector4(rays[ind][0],rays[ind][1],rays[ind][2]) 
        b = [a,da[0],da[1],da[2]]
        s = [0, rays[ind][3], rays[ind][4], rays[ind][5]]
        S.extend(s)
        m.extend(b)   
    return m, S

def polynomial_surface(point):
    x=point[0];y=point[1];z=point[2]
    dx,dy,dz = dvector4(x,y,z)
    xyz = vector4(x,y,z)
    H = ddvector4(x,y,z)
    S = [vector4(x,y,z),dx,dy,dz]
    c = [dot(v,vector4(x,y,z)),dot(v,dx),dot(v,dy),dot(v,dz)]
    for i in H:
        S.append(i)
        c.append(dot(v,i))
    A=[]
    b=[]
    for ind in idx:
        a=vector4(rays[ind][0],rays[ind][1],rays[ind][2])   
        da=dvector4(rays[ind][0],rays[ind][1],rays[ind][2]) 
        av = [a,da[0],da[1],da[2]]
        bv = [0, rays[ind][3], rays[ind][4], rays[ind][5]]
        A.extend(av)
        b.extend(bv)
    b = 2*dot(transpose(A),b)
    AtA = 2*dot(transpose(A),A)
    w = concatenate((b,c))
    m, _ = shape(AtA)
    p, q = shape(S)
    U = zeros((m+p,m+p))
    U[0:m,0:m]=AtA
    U[0:m,m:m+p]=transpose(S)
    U[m:m+p,0:m]=S
    x = solve(U,w)
    update_v(x[0:m])
    return

def standard_surface(point):
    m, S = matrix()
    A = dot(transpose(m),m)
    b = dot(transpose(m),S)
    v = solve(A, b)
    update_v(v)
    return

def update_surface(x,y,z):
    point = array([x,y,z])
    dif = point-center
    ratio = dot(dif,dif)/r**2
    if ratio > 0.25:
        idx = getSamples(point,r)
        #if len(idx)==0:
        #    return -1
        update_idx(idx)
        update_center(point)
        if ratio > 0.8: #update v
            standard_surface(point)
        else:
            polynomial_surface(point)
    return  0
    

def dfn4(x,y,z):
    update_surface(x,y,z)
##    print 'shape v',shape(v)
##    print 'shape dvector4(x,y,z)',shape(dvector4(x,y,z))
    return dot(dvector4(x,y,z),transpose(v))

@vectorize
def fn4(x,y,z):
    #if update_surface(x,y,z) == -1:
    #    return -1
    return dot(v,vector4(x,y,z))

def plot():
    ax = implicit.plot_implicit(fn4,bbox=(-5,5,-5,5,-5,5))
    Axes3D.scatter(ax,rays[:,0],rays[:,1],rays[:,2])
    plt.show()
