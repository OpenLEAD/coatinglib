from numpy import array, dot, cross, outer, eye, array_equal, sum, sqrt
from math import cos, sin, ceil, pi, isnan
from openravepy import IkFilterOptions
from abc import ABCMeta, abstractmethod

class KinBodyError(Exception):    
    def __init__(self):
        Exception.__init__(self, "object is not a KinBody.")

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
    p0 -- initial point.
    tol -- tolerance, stop criteria.
    """
    while True:
        df1 = s1.df(p0); f1 = s1.f(p0)
        df2 = s2.df(p0); f2 = s2.f(p0)

        crossdf = cross(df1,df2)
        if isnan(crossdf[0]) or isnan(crossdf[1]) or isnan(crossdf[2]):
            raise ValueError('Outside model')
        if array_equal(crossdf,array([0,0,0])):
            raise ValueError('Gradients are colinear')
        
        df1df2 = dot(df1,df2); df1df1 = dot(df1,df1); df2df2 = dot(df2,df2)
        beta = (-f1*df1df2+f2*df1df1)/(df1df2**2-df1df1*df2df2)
        alpha = (-f1*df2df2+f2*df1df2)/(df1df1*df2df2-df1df2*df1df2)
        dk = alpha*df1+beta*df2
        p1 = p0+dk
        if dot(dk,dk)<tol**2:
            norm = s1.df([p1[0],p1[1],p1[2]])
            norm /= sqrt(dot(norm,norm))
            p1 = array([p1[0],p1[1],p1[2],norm[0],norm[1],norm[2]])
            return p1
        else: p0=p1

def surfaces_tangent(ray, s2):
    """ Find the tangent of two surfaces in point/normal ray 

    Keyword arguments:
    ray -- 6x1 vector, point and normal in surface 1.
    s2 -- surface 2.
    """
    tan = cross(ray[3:6],s2.df(ray))
    tan = tan/sqrt(dot(tan,tan))   
    return tan

def hat(vec):
    """ Skew-symmetric matrix """
    hvec = array([[0, -vec[2], vec[1]],[vec[2],0,-vec[0]],[-vec[1], vec[0],0]])
    return hvec

def Rab(a,b):
    """ Matrix rotation from 'a' to 'b' """
    a = a/sqrt(dot(a,a))
    b = b/sqrt(dot(b,b))

    v = cross(a,b)
    sinab = sqrt(dot(v,v))
    vhat = hat(v)
    cosab = dot(a,b)

    if sinab == 0 and cosab == 1:
        R = eye(3)
    elif sinab == 0 and cosab == -1:
        R = -eye(3)
    else:    
        R = eye(3)+vhat+dot(vhat,vhat)*(1-cosab)/(sinab**2)
    return R

def Raxis(a,theta):
    """ Matrix rotation on 'a' axis, angle 'theta' """
    
    R = eye(3)*cos(theta) + hat(a)*sin(theta) + (1-cos(theta))*outer(a, a)
    return R

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
        p = array(p)
        return sum(p[:,0:3]*p[:,0:3],1)-self._Rn**2
        
    def find_iter(self, p0):
        self._Rn = sqrt(p0[0]**2+p0[1]**2+p0[2]**2)-self.coatingstep

    def criteria(self):
        return self._Rn>self.stopR

    def findnextparallel(self, p):
        d = sqrt(p[0]**2+p[1]**2+p[2]**2)
        n = ceil((d-self.stopR)/self.coatingstep)
        self._Rn = min(self.stopR+n*self.coatingstep, self._Rn)
        return

    def name(self):
        return 'sphere'

class Plane(IterSurface):
    """ Plane surface class. An object of this class can be
    a surface to be iterated and generate the coating trajectories.

    Keyword arguments:
    Rn0 -- height of the plane (z)
    stopR -- stop parameter of the iteration
    coatingstep -- iter step
    """

    def __init__(self, Rn0=3.770, stopR=1.59, coatingstep = 0.003):
        IterSurface.__init__(self, Rn0, stopR, coatingstep)

    def update(self):
        self._Rn = self._Rn-self.coatingstep

    def f(self, p):
        return p[2]-self._Rn

    def df(self, p):
        return array([0, 0, 1])

    def f_array(self, p):
        p = array(p)
        return p[:,2]-self._Rn

    def find_iter(self, p0):
        self._Rn = p0[2]-self.coatingstep

    def criteria(self):
        return self._Rn>self.stopR

    def findnextparallel(self, p):
        d = p[2]
        n = ceil((d-self.stopR)/self.coatingstep)
        self._Rn = min(self.stopR+n*self.coatingstep, self._Rn)
        return

    def name(self):
        return 'plane'
