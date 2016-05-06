from numpy import *
from scipy.spatial import KDTree
import time

b = load('blade_sampling_full/b_crop.npz')
r = load('blade_sampling_full/r_crop.npz')


b  = b['array']
r  = r['array']

approachrays = concatenate((b,r))

#approachrays = load('blade_sampling_full/t_crop.npz')
#approachrays = approachrays['array']
def tree_time(approachrays,r=0.05):
    T = KDTree(approachrays[:,0:3])
    rays = []
    N = len(approachrays)
    i=0
    t=time.time()
    initialTime = t
    I = ones((len(approachrays), 1), dtype=bool)
    while True:
        if I[i]:
            rays.append(approachrays[i])
            idx = T.query_ball_point(approachrays[i,0:3],r)
            idx = array(idx)
            idx = idx[idx>i]
            for j in idx:I[j]=False
        if i%10000==0:
            print str(i)+'/'+str(N)
            print 'elapsed time= '+str(time.time()-t)
            t=time.time()
        i+=1
        if i==len(approachrays):break
    print 'total Time = '+str(time.time()- initialTime)
    return rays
approachrays = tree_time(approachrays,r=0.05)
savez_compressed('blade_sampling_full/blade_crop_fast.npz', array=approachrays)
