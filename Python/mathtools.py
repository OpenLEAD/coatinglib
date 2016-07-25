from numpy import array, dot
from math import cos, sin

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
