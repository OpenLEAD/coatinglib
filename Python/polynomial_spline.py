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
reachableRays  = reachableRays['array']
AllreachableRays = load('blade_sampling/blade_crop_fast2.npz')
AllreachableRays  = AllreachableRays['array']
rays = concatenate((reachableRays,AllreachableRays))


v=[]
r=0.4
center=[500,0,0]
idx=[]
handles2=[]
env2=None
index=-1
# ======= PESO
s = 0.11
c1 = sqrt(2*pi*s**2)
delta = 2*s**2
def evaluate_gauss(point,mean):
    dif = point-mean
    k = -dot(dif,dif)/delta
    return exp(k)/c1

def compute_W(point,rays,idx):
    w = []
    for i in idx:
        w.append(evaluate_gauss(rays[i][0:3],point))
    return w
# =======

def set_handles(handles,env):
    global handles2, env2
    handles2 = handles
    env2 = env
    return

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

def dv_vector4(x,y,z):
    x2=x**2;y2=y**2;z2=z*z;xy=x*y
    x3=x2*x;y3=y2*y;z3=z2*z;yz=y*z
    x4=x3*x;y4=y3*y;z4=z3*z;xz=x*z
    
    
    a = [1, z, z2, z3, z4, y, yz, y*z2, y*z3, y2, y2*z, y2*z2, y3, y3*z, y4, x, xz, x*z2, x*z3, xy, xy*z, #20
            xy*z2, x*y2, xy*yz, x*y3, x2, x2*z, x2*z2, x2*y, x2*yz, x2*y2, x3, x3*z, x3*y, x4]#34
    
    dx = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, z, z2, z3, y, yz, a[7], y2, a[10], y3,
          2*x, 2*xz, 2*a[17], 2*xy, 2*a[20], 2*a[22], 3*x2, 3*a[26], 3*a[28], 4*x3]
    dy = [0, 0, 0, 0, 0, 1, z, z2, z3, 2*y, 2*yz, 2*a[7], 3*y2, 3*a[10], 4*y3, 0, 0,
          0, 0, x, xz, a[17], 2*xy, 2*a[20], 3*a[22], 0, 0, 0, x2, a[26], 2*a[28], 0, 0, x3, 0]
    dz = [0, 1, 2*z, 3*z2, 4*z3, 0, y, 2*yz, 3*a[7], 0, y2, 2*a[10], 0, y3, 0, 0, x, 2*xz, 3*a[17], 0,
          xy, 2*a[20], 0, a[22], 0, 0, x2, 2*a[26], 0, a[28], 0, 0, x3, 0, 0]
    da=[dx,dy,dz]
    return [a,da]

def ddvector4(x,y,z):
    x2=x*x;y2=y*y;z2=z*z
    xy=x*y;yz=y*z;xz=x*z
    
    dxdx = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 2, 2*z, 2*z2, 2*y, 2*yz, 2*y2, 6*x, 6*xz, 6*xy, 12*x2]
    dxdy = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, z, z2,
            2*y, 2*yz, 3*y2, 0, 0, 0, 2*x, 2*xz, 4*xy, 0, 0, 3*x2, 0]
    dxdz = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2*z, 3*z2, 0,
            y, 2*yz, 0, y2, 0, 0, 2*x, 4*xz, 0, 2*xy, 0, 0, 3*x2, 0, 0]
    #dydx = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, z, z**2, 2*y, 2*y*z, 3*y**2, 0, 0, 0, 2*x, 2*x*z, 4*x*y, 0, 0, 3*x**2, 0]
    dydy = [0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2*z, 2*z2, 6*y, 6*yz, 12*y2, 0, 0,
            0, 0, 0, 0, 0, 2*x, 2*xz, 6*xy, 0, 0, 0, 0, 0, 2*x2, 0, 0, 0, 0]
    dydz = [0, 0, 0, 0, 0, 0, 1, 2*z, 3*z2, 0, 2*y, 4*yz, 0, 3*y2, 0, 0, 0, 0, 0,
            0, x, 2*xz, 0, 2*xy, 0, 0, 0, 0, 0, x2, 0, 0, 0, 0, 0]
    #dzdx = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2*z, 3*z**2, 0, y, 2*y*z, 0, y**2, 0, 0, 2*x, 4*x*z, 0, 2*x*y, 0, 0, 3*x**2, 0, 0]
    #dzdy = [0, 0, 0, 0, 0, 0, 1, 2*z, 3*z**2, 0, 2*y, 4*y*z, 0, 3*y**2, 0, 0, 0, 0, 0, 0, x, 2*x*z, 0, 2*x*y, 0, 0, 0, 0, 0, x**2, 0, 0, 0, 0, 0]
    dzdz = [0, 0, 2, 6*z, 12*z2, 0, 0, 2*y, 6*yz, 0, 0, 2*y2, 0, 0, 0, 0, 0,
            2*x, 6*xz, 0, 0, 2*xy, 0, 0, 0, 0, 0, 2*x2, 0, 0, 0, 0, 0, 0, 0]
    return [dxdx,dxdy,dxdz,dydy,dydz,dzdz]

def make_matrix():
    m = []
    S = []
    for ind in idx:
        [a,da]=dv_vector4(rays[ind][0],rays[ind][1],rays[ind][2])   
        b = [a,da[0],da[1],da[2]]
        s = [0, rays[ind][3], rays[ind][4], rays[ind][5]]
        S.extend(s)
        m.extend(b)   
    return m, S

def polynomial_surface(point):
    [a,da]=dv_vector4(point[0],point[1],point[2]) 
    H = ddvector4(point[0],point[1],point[2])
    S = [a,da[0],da[1],da[2]]
    c = [dot(v,a),dot(v,da[0]),dot(v,da[1]),dot(v,da[2])]
    for i in H:
        S.append(i)
        c.append(dot(v,i))
    A=[]
    b=[]
    for ind in idx:
        [a,da]=dv_vector4(rays[ind][0],rays[ind][1],rays[ind][2])
        a=vector4(rays[ind][0],rays[ind][1],rays[ind][2])   
        av = [a,da[0],da[1],da[2]]
        bv = [0, rays[ind][3], rays[ind][4], rays[ind][5]]
        A.extend(av)
        b.extend(bv)
    At = transpose(A)
    
    # ===============================
    # Polynomial Spline with weights
    w = compute_W(point,rays,idx)
    A=array(A)
    for i in range(0,len(A)):
        A[i]*=w[i/4]
    # ===============================
    
    b = 2*dot(transpose(A),b)
    AtA = 2*dot(At,A)
    w = concatenate((b,c))
    m, _ = shape(AtA)
    p, _ = shape(S)
    U = zeros((m+p,m+p))
    U[0:m,0:m]=AtA
    U[0:m,m:m+p]=transpose(S)
    U[m:m+p,0:m]=S
    x = solve(U,w)
    update_v(x[0:m])
    return

def standard_surface(point):
    m, S = make_matrix()
    A = dot(transpose(m),m)
    b = dot(transpose(m),S)
    v = solve(A, b)
    update_v(v)
    return

def update_surface(x,y,z):
    global handles, index
    point = array([x,y,z])
    dif = point-center
    ratio = dot(dif,dif)/r**2
    if ratio > 0.25: #0.01
        idx = getSamples(point,r)
        #if len(idx)==0:
        #    return -1
        update_idx(idx)
        update_center(point)
        if ratio > 0.85: #update v #0.25
            print 'arroe! (x,y,z) - ', point
            standard_surface(point)
        else:
            print 'recalculando...'
            #if index>-1: handles2.pop(index-1) 
            #handles2.append(env2.plot3(points=rays[idx,0:3],pointsize=5,colors=array((0,0,0))))
            #index=len(handles2)
            polynomial_surface(point)

        maxi=0
        ponto=[]
        for ind in idx:
            tempfn4 = abs(dot(v,vector4(rays[ind][0],rays[ind][1],rays[ind][2])))
            if maxi<=tempfn4:
                maxi=tempfn4
                ponto = rays[ind]
        #print 'Pior ponto =', ponto[0:3], ' fn4(x,y,z) = ',maxi
        #handles2.append(env2.plot3(points=ponto[0:3],pointsize=5,colors=array((1,0,0))))
        #wait = input("PRESS ENTER TO CONTINUE.")    
                
        #print v    
    return  0
    

def dfn4(x,y,z):
    update_surface(x,y,z)
    return dot(dvector4(x,y,z),transpose(v))

@vectorize
def fn4(x,y,z):
    update_surface(x,y,z)
    return dot(v,vector4(x,y,z))

def plot():
    ax = implicit.plot_implicit(fn4,bbox=(-5,5,-5,5,-5,5))
    Axes3D.scatter(ax,rays[:,0],rays[:,1],rays[:,2])
    plt.show()

def optmizeTan(p0, pnew, tol):
    def func(P):
        a=P-pnew
        return dot(a,a)
    def func_deriv(P):
        return 2*(P-pnew)
    def consfunc(P):
        return fn4(P[0],P[1],P[2])
    def consfunc_deriv(P):
        return dfn4(P[0],P[1],P[2])
    cons = ({'type':'eq',
             'fun': consfunc,
            'jac':consfunc_deriv})
    res = minimize(func, p0,jac=func_deriv,constraints=cons, method='SLSQP',tol=tol, options={'disp': False})
    return res    
