from sympy import *

x = Symbol('x')
y = Symbol('y')
z = Symbol('z')

def derivate_vector(vector,variable):
    dvector=[]
    for i in vector:
        if type(i)==int:
            dvector.append(0)
        else:
            dvector.append(i.diff(variable))
    return dvector		
