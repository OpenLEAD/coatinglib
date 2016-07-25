from numpy import array, dot
from math import cos, sin, sqrt

def Rotate(obj, theta, axis):
    """ Rotate object in x, y, or z axis (world frame)

    Keyword arguments:
    obj -- object.
    theta -- angle in rad.
    axis -- string x, y, or z.
    """
    
    if axis=='z':
        return RzM(obj,theta)
    elif axis=='x':
        return RxM(obj,theta)
    elif axis=='y':
        return RyM(obj,theta)
    elif axis=='':
        return obj
    else:
        print "mathtools::Rotate wrong axis"
        return obj

def T(theta, axis):
    """ Transformation matrix axis

    Keyword arguments:
    theta -- angle in rad.
    axis -- string x, y, or z.
    """
    
    if axis=='z':
        return Tz(theta)
    elif axis=='x':
        return Tx(theta)
    elif axis=='y':
        return Ty(theta)
    elif axis=='':
        return eye(4)
    else:
        print "mathtools::T wrong axis"
        return eye(4)    
        
def Tx(theta):
    """ Transformation matrix x

    theta -- angle in rad.
    """

    T = array([[1, 0, 0, 0],
                   [0, cos(theta), -sin(theta), 0],
                   [0, sin(theta), cos(theta), 0],
                   [0, 0, 0, 1]])
    return T

def Ty(theta):
    """ Transformation matrix y

    theta -- angle in rad.
    """

    T = array([[cos(theta), 0, sin(theta), 0],
                   [0, 1, 0, 0],
                   [-sin(theta), 0, cos(theta), 0],
                   [0, 0, 0, 1]])
    return T

def Tz(theta):
    """ Transformation matrix z

    theta -- angle in rad.
    """

    T = array([[cos(theta), -sin(theta), 0, 0],
                   [sin(theta), cos(theta), 0, 0],
                   [0, 0, 1, 0],
                   [0, 0, 0, 1]])
    return T
    
def RzM(obj,theta):
    """ Rotate object in Z axis (world frame)

    Keyword arguments:
    obj -- object.
    theta -- angle in rad.
    """

    T = dot(array([[cos(theta), -sin(theta), 0, 0],
                   [sin(theta), cos(theta), 0, 0],
                   [0, 0, 1, 0],
                   [0, 0, 0, 1]]),
            obj.GetTransform())
    
    obj.SetTransform(T)
    return obj

def RxM(obj,theta):
    """ Rotate object in X axis (world frame)

    Keyword arguments:
    obj -- object.
    theta -- angle in rad.
    """
    T = dot(array([[1, 0, 0, 0],
                   [0, cos(theta), -sin(theta), 0],
                   [0, sin(theta), cos(theta), 0],
                   [0, 0, 0, 1]]),
            obj.GetTransform())
    
    obj.SetTransform(T)
    return obj

def RyM(obj,theta):
    """ Rotate object in Y axis (world frame)

    Keyword arguments:
    obj -- object.
    theta -- angle in rad.
    """
    T = dot(array([[cos(theta), 0, sin(theta), 0],
                   [0, 1, 0, 0],
                   [-sin(theta), 0, cos(theta), 0],
                   [0, 0, 0, 1]]),
            obj.GetTransform())
    
    obj.SetTransform(T)
    return obj

def curvepoint(s1, s2, p0, tol=1e-4):
    """ Find a point near p0 which is in f_1 and f_2 (intersection between surfaces) 

    Keyword arguments:
    s1 -- surface 1.
    s2 -- surface 2.
    p0 -- initial point.
    tol -- tolerance, stop criteria.
    """
    
    while True:
        df1 = s1.df(p0); f1 = s1.f(p0)
        df2 = s2.df(p0); f2 = s2.f(p0)
        df1df2 = dot(df1,df2); df1df1 = dot(df1,df1); df2df2 = dot(df2,df2)
        beta = (-f1+f2*df1df1/df1df2)/(df1df2-df1df1*df2df2/df1df2)
        alpha = (-f1+f2*df1df2/df2df2)/(df1df1-df1df2*df1df2/df2df2)
        dk = alpha*df1+beta*df2
        p1 = p0+dk
        if sqrt(dot(dk,dk))<tol:
            return p1
        else: p0=p1

class sphere:
    """ Sphere surface class. An object of this class can be
    a surface to be iterated and generate the coating trajectories.

    Keyword arguments:
    Rn0 - first radius of the sphere
    """

    def __init__(self, Rn0):
        self._Rn = Rn0

    def update(self):
        self._Rn = self._Rn-0.003

    def f(self, p):
        return p[0]**2+p[1]**2+p[2]**2-self._Rn**2

    def df(self, p):
        return array([2*p[0], 2*p[1], 2*p[2]])

    def find_iter(self, p0):
        self._Rn = sqrt(p0[0]**2+p0[1]**2+p0[2]**2)-0.003
