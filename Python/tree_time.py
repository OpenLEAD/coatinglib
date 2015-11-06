from numpy import *
from scipy.spatial import KDTree
import time


approachrays = load('blade_faro2.npz')
approachrays = approachrays['array']

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
        idx = T.query_ball_point(approachrays[i,0:3],r=0.01)
        idx = array(idx)
        idx = idx[idx>i]
        for j in idx:I[j]=False
    if i%1000==0:
        print str(i)+'/'+str(N)
        print 'elapsed time= '+str(time.time()-t)
        t=time.time()
    i+=1
    if i==len(approachrays):break
savez_compressed('blade_faro_fast2.npz', array=rays)
print 'total Time = '+str(time.time()- initialTime)
