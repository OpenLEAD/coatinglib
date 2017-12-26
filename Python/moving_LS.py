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


Rays1 = load('coated_points/pos4_black.npz')
Rays1  = Rays1['array']
Rays2 = load('coated_points/pos4_blue.npz')
Rays2  = Rays2['array']
rR = array(concatenate((Rays1,Rays2)))

## GAUSS PARAMETERS
s = 5
c1 = sqrt(2*pi*s**2)
delta = 2*s**2
r=s*sqrt(10*log(10))

##u2 = (rR[1000]+rR[1001])/2
##u2[3:6]/=sqrt(dot(u2[3:6],u2[3:6]))

# TREE TIME
Tree = KDTree(rR[:,0:3])

def update_normal_const(sigma):
    global s, c1, delta, r
    s = sigma
    c1 = sqrt(2*pi*s**2)
    delta = 2*s**2
    r=s*sqrt(10*log(10))
    return


def vector4(x,y,z):
    x2=x**2;y2=y**2;z2=z*z;xy=x*y
    x3=x2*x;y3=y2*y;z3=z2*z;yz=y*z
    x4=x3*x;y4=y3*y;z4=z3*z
    return [1, z, z2, z3, z4, y, y*z, y*z2, y*z3, y2, y2*z, y2*z2, y3, y3*z, y4, x, x*z, x*z2, x*z3, xy, xy*z,
            xy*z2, x*y2, xy*yz, x*y3, x2, x2*z, x2*z2, x2*y, x2*yz, x2*y2, x3, x3*z, x3*y, x4]

def dvector4(x,y,z):
    dx = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, z, z**2, z**3, y, y*z, y*z**2, y**2, y**2*z, y**3, 2*x, 2*x*z, 2*x*z**2, 2*x*y, 2*x*y*z, 2*x*y**2, 3*x**2, 3*x**2*z, 3*x**2*y, 4*x**3]
    dy = [0, 0, 0, 0, 0, 1, z, z**2, z**3, 2*y, 2*y*z, 2*y*z**2, 3*y**2, 3*y**2*z, 4*y**3, 0, 0, 0, 0, x, x*z, x*z**2, 2*x*y, 2*x*y*z, 3*x*y**2, 0, 0, 0, x**2, x**2*z, 2*x**2*y, 0, 0, x**3, 0]
    dz = [0, 1, 2*z, 3*z**2, 4*z**3, 0, y, 2*y*z, 3*y*z**2, 0, y**2, 2*y**2*z, 0, y**3, 0, 0, x, 2*x*z, 3*x*z**2, 0, x*y, 2*x*y*z, 0, x*y**2, 0, 0, x**2, 2*x**2*z, 0, x**2*y, 0, 0, x**3, 0, 0]
    a=[dx,dy,dz]
    return a

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

def evaluate_gauss(point,mean):
    dif = point-mean
    k = -dot(dif,dif)/delta
    return exp(k)/c1

def compute_W(point,rays):
    w = []
    for i in range(0,len(rays)):
        w.append(evaluate_gauss(rays[i][0:3],point))
    return w

def multiply_by_w(w,m):
    for i in range(0,len(m)):
        m[i]*=w[i/4]
    return wm    

def compute_dv(point,rays,v,AtwA,A,m,S):
 #   print 'len(m)=',len(m)
  #  t = time.time()
    #W = rays
    Wxh = []
    Wyh = []
    Wzh = []
    for i in range(0,len(m)):
        h = m[i]*2.0/delta
        #W = rays[i,:]*h
        Wxh.append(rays[i/4][0]*h)
        Wyh.append(rays[i/4][1]*h)
        Wzh.append(rays[i/4][2]*h)
  #  print 'for_computedv_time = '+str(time.time()- t)    
  #  t = time.time()    
    dif = S-dot(A,v)

    dvx = solve(AtwA,dot(transpose(Wxh),dif))
    dvy = solve(AtwA,dot(transpose(Wyh),dif))
    dvz = solve(AtwA,dot(transpose(Wzh),dif))
    a=[dvx,dvy,dvz]
  #  print 'solve_computedv_time = '+str(time.time()- t)
    return a
        
def matrix(rays):
    m = []
    S = []
    for ind in range(0,len(rays)):
##        t = time.time()
 #       a=vector4(rays[ind][0],rays[ind][1],rays[ind][2])
##        print 'vector4_time = '+str(time.time()- t)
##        t=time.time()
        
 #       da=dvector4(rays[ind][0],rays[ind][1],rays[ind][2])
##        print 'dvector4_time = '+str(time.time()- t)
##        t=time.time()
        [a,da] = dv_vector4(rays[ind][0],rays[ind][1],rays[ind][2])
        b = [a,da[0],da[1],da[2]]
        s = [0, rays[ind][3], rays[ind][4], rays[ind][5]]
        S.extend(s)
        m.extend(b)
##        print 'extend_time = '+str(time.time()- t)
##        t=time.time()
        
    return m, S

def polynomial_surface(point,rays):
    A, S = matrix(rays)
    w = compute_W(point,rays)
     

    S = array(S)
    A = array(A)
    m = copy.copy(A)

    At=transpose(m)
    for i in range(0,len(m)):
        m[i]*=w[i/4] # m = W*A   
    
    AtwA = dot(At,m)
    b = dot(transpose(m),S)
  
    v = solve(AtwA, b)
    return v, AtwA, A, m, S       

def dpolynomial(point,rays):
##    t = time.time()
    idx = Tree.query_ball_point(point,r)
    v, AtwA, A, m, S  = polynomial_surface(point,rays)
    dv = compute_dv(point,rays,v,AtwA,A,m,S)
    B = vector4(point[0],point[1],point[2])
    dB = dvector4(point[0],point[1],point[2])
    dP = dot(dv,B)+dot(dB,v)
##    print 'dpolynomial_time = '+str(time.time()- t)
    return dP

@vectorize
def fn4(x,y,z):
    point=array([x,y,z])
    v, _, _, _, _ = polynomial_surface(point,rR)
    return dot(v,vector4(x,y,z))

def fn4_1(P):
    x,y,z=P
    return [fn4(x,y,z),0,0]

def dfn4(x,y,z):
    point = [x,y,z]
    return dpolynomial(point,rR)

def plot():
    ax = implicit.plot_implicit_save(fn4)
    #Axes3D.scatter(ax,rR[:,0],rR[:,1],rR[:,2])
    #plt.show()

def optmizeTan(p0, pnew, tol):
    def func(P):
        a=P-pnew
        return dot(a,a)
    def func_deriv(P):
        return 2*(P-pnew)
    def consfunc(P):
        return fn4(P[0],P[1],P[2])
    def consfunc_deriv(P):
        return dpolynomial(P,rR)
    cons = ({'type':'eq',
             'fun': consfunc,
            'jac':consfunc_deriv})
    res = minimize(func, p0,jac=func_deriv,constraints=cons, method='SLSQP',tol=tol, options={'disp': False})
    return res
