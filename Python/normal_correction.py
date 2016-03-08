from scipy.spatial import KDTree
from numpy import *
import time

approachrays = load('blade_sampling/blade_pcl.npz')
approachrays = approachrays['array']

def normal_correction(approachrays,r=0.05):
    T = KDTree(approachrays[:,0:3])
    rays = []
    N = len(approachrays)
    t=time.time()
    initialTime = t
    for i in range(0,len(approachrays)):
        rays.append(approachrays[i])
        idx = T.query_ball_point(approachrays[i,0:3],r)
        for ind in idx:
            if dot(approachrays[i,3:6],approachrays[ind,3:6])<0:
                approachrays[ind,3:6]=-approachrays[ind,3:6]
        if i%10000==0:
            print str(i)+'/'+str(N)
            print 'elapsed time= '+str(time.time()-t)
            t=time.time()
    return rays

def main():
    rays = normal_correction(approachrays,r=0.05)
    return rays
