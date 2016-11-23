from numpy import array, dot, cross, outer, eye, sum, sqrt, ones
from numpy import random, transpose, zeros, linalg, multiply, linspace, power
from numpy import ndindex, linspace, power,  ceil, floor, einsum
from numpy import cos as npcos
from numpy import sin as npsin
from math import cos, sin, pi, isnan
from abc import ABCMeta, abstractmethod
from copy import copy
from collections import deque

class KinBodyError(Exception):    
    def __init__(self):
        Exception.__init__(self, "object is not a KinBody.")

def annulus_distribution(N,r,R, origin = None, dim = 2):
    """
    return N samples from a uniform distribution over a
    annulus with inner radius r and outer radius R.

    """
    if origin is None:
        origin = zeros((1,dim))
    
    direction = random.randn(N,dim)
        
    rho = random.rand(N)
    rho = power((R**dim - r**dim)*rho + r**dim,1./dim)
    
    samples = einsum('ij,i->ij',direction,rho/linalg.norm(direction, axis=1))

    return samples + origin


def fast_poisson_disk(N, r, limits, k = 30):
    n = len(limits)
    cellfreq = sqrt(n)/r
    delta = dot(limits,[-1, 1])
    gridsize = ceil(delta*cellfreq)
    grid = - ones(gridsize.astype('int')).astype('int') 
    x0 = random.rand(n)
    grid[tuple(floor(x0*gridsize).astype('int'))] = 0
    x = [limits[:,0]+delta*x0]
    activelist = [0]
    counter = 0
    while len(activelist)>0:
        i = activelist[random.randint(len(activelist))]
        samples = annulus_distribution(k,r,2*r,x[i],n)
        placed = False
        for sample in samples:
            grid_place = floor((sample-limits[:,0])*gridsize/delta).astype('int')
            if any(grid_place < 0) or any(grid_place >= grid.shape):
                continue

            for idx in ndindex(tuple([3]*n)):
                try:
                    point_index = grid[tuple(array(idx)-1+grid_place)]
                except IndexError:
                    continue
                if point_index == -1:
                    continue
                if linalg.norm(sample - x[point_index]) < r:
                    break
            else:
                grid[tuple(grid_place)] = len(x)
                activelist += [len(x)]
                x += [sample]
                placed = True
        if not placed:
            del activelist[activelist.index(i)]
    return x
            

def central_difference(turbine, joints_trajectory, trajectory_index):

    # 'j' for joints, so it doesnt clumpsy the equations
    j = deque(joints_trajectory)
    j.rotate(-trajectory_index)

    h = turbine.config.model.trajectory_step
    
    # Joints velocity - Central Difference (h**6 order error)
    w = ( (j[3]-j[-3]) + 9*(j[-2]-j[2]) + 45*(j[1]-j[-1]) )/(60.0*h)

    # Joints acceleration - Central Difference (h**6 order error)
    alpha = ( 2*(j[-3]+j[3]) - 27*(j[-2]+j[2]) + 270*(j[-1]+j[1]) - 490*j[0] )/(180.0*h**2)

    return w, alpha

def curvepoint(s1, s2, p0, tol=1e-4):
    """
    Find a point on s1 and s2 (intersection between surfaces) near p0.
    The algorithm is described in the paper:
    'Tracing surface intersections - C.L. Bajaj, C.M. Hoffmann, R.E. Lynch,
    and J.E.H. Hopcroft.'
    Available at:
    https://www.cs.purdue.edu/homes/cmh/distribution/papers/Geometry/geo6.pdf

    Keyword arguments:
    s1 -- surface 1.
    s2 -- surface 2.
    p0 -- initial point (3D - without normal vector).
    tol -- tolerance, stop criteria.
    """
    
    tol2 = tol**2
    while True:
        df1 = s1.df(p0); f1 = s1.f(p0)
        df2 = s2.df(p0); f2 = s2.f(p0)

        crossdf = cross(df1,df2)
        if isnan(crossdf[0]) or isnan(crossdf[1]) or isnan(crossdf[2]):
            raise ValueError('Outside model')
        if (abs(crossdf - array([0,0,0]))<=1e-5).all():
            raise ValueError('Gradients are colinear')
        
        df1df2 = dot(df1,df2); df1df1 = dot(df1,df1); df2df2 = dot(df2,df2)
        beta = (-f1*df1df2+f2*df1df1)/(df1df2**2-df1df1*df2df2)
        alpha = (-f1*df2df2+f2*df1df2)/(df1df1*df2df2-df1df2*df1df2)
        dk = alpha*df1+beta*df2
        p1 = p0+dk
        if dot(dk,dk)<tol2:
            grad = s1.df(p1)
            grad = grad/sqrt(dot(grad,grad))
            p1 = array([p1[0],p1[1],p1[2],
                        grad[0],grad[1],grad[2]])
            return p1
        else: p0=p1

def surfaces_tangent(ray, s2):
    """ Find the tangent of two surfaces in point/normal ray 

    Keyword arguments:
    ray -- 6x1 vector, point and normal in surface 1.
    s2 -- surface 2.
    """
    
    tan = cross(ray[3:6],s2.df(ray))
    tan = tan/linalg.norm(tan)
    return tan

def hat(vec):
    """ Skew-symmetric matrix """
    
    hvec = array([[0, -vec[2], vec[1]],
                  [vec[2],0,-vec[0]],
                  [-vec[1], vec[0],0]])
    return hvec

def Rab(a,b):
    """ Matrix rotation from 'a' to 'b' """
    a = a/linalg.norm(a)
    b = b/linalg.norm(b)

    v = cross(a,b)
    sinab = sqrt(dot(v,v))
    vhat = hat(v)
    cosab = dot(a,b)

    if cosab == -1:
        return Raxis(compute_perpendicular_vector(a), pi)
    else:    
        R = eye(3)+vhat+dot(vhat,vhat)*1.0/(1+cosab)
    return R

def Raxis(a,theta):
    """ Matrix rotation on 'a' axis, angle 'theta' """
    
    a = a/linalg.norm(a)    
    R = eye(3)*cos(theta) + hat(a)*sin(theta) + (1-cos(theta))*outer(a, a)
    return R

def rotate_trajectories(turbine, trajectories, T=[]):
    """
    Rotate all points in trajectories and return the trajectories.

    Keyword arguments:
    T -- the homogeneous transform matrix 4x4.
    """
    
    if len(T)==0:
        T = turbine.blades[0].GetTransform()

    R = T[0:3,0:3]
    for i in range(0,len(trajectories)):
        traj = array(trajectories[i])
        Ra = zeros(traj.shape)
        
        Ra[:,0:3] = dot(traj[:,0:3], transpose(R))
        Ra[:,3:6] = dot(traj[:,3:6], transpose(R))
        
        trajectories[i] = Ra.tolist()
    return trajectories
        
    

def compute_perpendicular_vector(vector_1):
    """
    Compute a perpendicular vector of a given vector.

    Keyword arguments:
    vector_1 -- given vector.
    """
    
    vector_1 = vector_1/linalg.norm(vector_1)
    vector_2 = copy(vector_1)
    while abs(dot(vector_1,vector_2)-1)<=1e-6:
        vector_2 = vector_1+random.uniform(1,-1,3)
        vector_2 = vector_2/linalg.norm(vector_2)                  
    return cross(vector_1,vector_2)

def distance_point_line_3d(x1, x2, x0):
    """
    A line in three dimensions specified by two points x1 and x2.
    The function calculates the squared distance between a point
    x0 and the line.
    
    keyword arguments:
    x1 -- (x,y,z) a point of the line
    x2 -- (x,y,z) a second point of the line
    x0 -- point not in the line to compute distance
    """

    return linalg.norm(cross(x0-x1,x0-x2))/linalg.norm(x2-x1)

    
class IterSurface:
    """ Inheritable class to surfaces that can be iterated and generate the coating
    trajectories.

    Keyword arguments:
    Rn0 -- initial parameter of the iteration
    stopR -- stop parameter of the iteration
    coatingstep -- iter step
    """

    __metaclass__ = ABCMeta
    
    def __init__(self, Rn0=3.770, stopR=1.59, coatingstep = 0.003):
        self._Rn0 = Rn0
        self._Rn = Rn0
        self.stopR = stopR
        self.coatingstep = coatingstep

    def update(self):
        self._Rn = self._Rn-self.coatingstep

    @abstractmethod
    def f(self, p):
        pass 

    @abstractmethod
    def df(self, p):
        pass

    @abstractmethod
    def f_array(self, p):
        pass 

    @abstractmethod
    def find_iter(self, p0):
        pass

    @abstractmethod
    def criteria(self):
        pass 

    @abstractmethod
    def findnextparallel(self, p):
        pass

    @abstractmethod
    def name(self, p):
        pass


class Sphere(IterSurface):
    """ Sphere surface class. An object of this class can be
    a surface to be iterated and generate the coating trajectories.

    Keyword arguments:
    Rn0 -- initial parameter of the sphere
    stopR -- stop parameter of the iteration
    coatingstep -- iter step
    """

    def __init__(self, Rn0=3.770, stopR=1.59, coatingstep = 0.003):
        IterSurface.__init__(self, Rn0, stopR, coatingstep)
        
    def f(self, p):
        return p[0]**2+p[1]**2+p[2]**2-self._Rn**2

    def df(self, p):
        return array([2*p[0], 2*p[1], 2*p[2]])

    def f_array(self, p):
        if len(p)==1: return self.f(p[0])
        p = array(p)
        return sum(p[:,0:3]*p[:,0:3],1)-self._Rn**2
        
    def find_iter(self, p0):
        self._Rn = sqrt(p0[0]**2+p0[1]**2+p0[2]**2)

    def criteria(self):
        return self._Rn>self.stopR

    def findnextparallel(self, p):
        d = sqrt(p[0]**2+p[1]**2+p[2]**2)
        n = ceil((d-self.stopR)/self.coatingstep)
        self._Rn = min(self.stopR+n*self.coatingstep, self._Rn)
        return

    def name(self):
        return 'Sphere'

class Plane(IterSurface):
    """ Plane surface class. An object of this class can be
    a surface to be iterated and generate the coating trajectories.

    Keyword arguments:
    Rn0 -- height of the plane (z)
    stopR -- stop parameter of the iteration
    coatingstep -- iter step
    """

    def __init__(self, Rn0=1, stopR=-4, coatingstep = 0.003):
        IterSurface.__init__(self, Rn0, stopR, coatingstep)

    def f(self, p):
        return p[2]-self._Rn

    def df(self, p):
        return array([0, 0, 1])

    def f_array(self, p):
        if len(p)==1: return self.f(p[0])
        p = array(p)
        return p[:,2]-self._Rn

    def find_iter(self, p0):
        self._Rn = p0[2]

    def criteria(self):
        return self._Rn>self.stopR

    def findnextparallel(self, p):
        d = p[2]
        n = ceil((d-self.stopR)/self.coatingstep)
        self._Rn = min(self.stopR+n*self.coatingstep, self._Rn)
        return

    def name(self):
        return 'Plane'

_backdif2 = array([[1.0, -1.0],
                   [0,    0  ]])

_backdif3 = array([[3.0/2,    -2.0,   1.0/2 ],
                   [1.0,      -2.0,   1.0   ]])

_backdif4 = array([[11.0/6,   -3.0,   3.0/2,  -1.0/3],
                   [2.0,      -5.0,   4.0,    -1.0  ]])

_backdif5 = array([[25.0/12,  -4.0,    3.0,      -4.0/3,  1.0/4],
                   [35.0/12,  -26.0/3, 19.0/2,   -14.0/3, 11.0/12]])

_backdif6 = array([[137.0/60, -5.0,     5.0,      -10.0/3, 5.0/4,   -1.0/5],
                   [15.0/4,   -77.0/6,  107.0/6,  -13.0,   61.0/12, -5/6]])

_backdif7 = array([[49.0/20,  -6.0,     15.0/2,  -20.0/3,  15.0/4, -6.0/5,  1.0/6],
                   [203.0/45, -87.0/5,	117.0/4, -254.0/9, 33.0/2, -27.0/5, 137.0/180]])

_backdiflist = (_backdif2,
                _backdif3,
                _backdif4,
                _backdif5,
                _backdif6,
                _backdif7)

def backward_difference(turbine, joints_trajectory):
    """
    joints_trajectory - joints_trajectory[0] is the most recent
    """
    joints = array(joints_trajectory[:7])
    h = turbine.config.model.trajectory_step
    size = len(joints)

    w, alpha = h**(size-size/2) * dot(_backdiflist[size-2],
                                      h**(size/2) * joints)
    
    return w, alpha

def partial_backward_difference(turbine, old_joints_trajectory):
    """
    old_joints_trajectory - old_joints_trajectory[0] is the most recent past value
    Compute from old_joints_trajectory the partial backward diffeference, before adding the last point
    """
    joints = array(old_joints_trajectory[:6])
    h = turbine.config.model.trajectory_step
    size = len(joints)+1

    w_part, alpha_part = dot(_backdiflist[size-2][:,1:],
                                      h**(size/2) * joints)
    
    return w_part, alpha_part, size

def update_backward_difference(turbine, actual_joint, w_part, alpha_part, size):
    """
    joints_trajectory - joints_trajectory[0] is the most recent
    """
    h = turbine.config.model.trajectory_step

    w_alpha = h**(size-size/2)* (transpose((w_part, alpha_part)) +
                                  multiply.outer(h**(size/2) * actual_joint,
                                                 _backdiflist[size-2][:,0]))
        

    return w_alpha[:,:,0],w_alpha[:,:,1]
