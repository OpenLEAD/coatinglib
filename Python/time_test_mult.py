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

a=range(0,1000000)
b=range(0,1000000)
c=ones(1000000)*2
t=time.time()
for i in range(0,len(a)):
    b[i]=c[i]+a[i]
print 'sum_time = '+str(time.time()- t)

t=time.time()
for i in range(0,len(a)):
    b[i]=2.0*a[i]
print 'mult_time = '+str(time.time()- t)
