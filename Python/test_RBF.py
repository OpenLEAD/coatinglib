import RBF_c
from numpy import load, array
from os import remove

try:
    points = load('blade_sampling_full/blade_crop_ufast.npz')
    points  = points['array']
    test_point = points[0][0:3]+array([0.1,0.1,0.1])
except:
    raise ValueError('Points are not in the right folder')

try:
    remove('./RBF/jiraublade_r3_points.npz')
    remove('./RBF/jiraublade_r3_w.npz')
except:
    None

try:
    rbf = RBF_c.RBF('jiraublade', 'r3', points)
except: raise

try: rbf.make()
except: raise

try: rbf.f(test_point)
except: raise

try: rbf.df(test_point)
except: raise

try:
    del rbf
    rbf = RBF_c.RBF('jiraublade', 'r3')
except: raise    

print "Test succeeded"
