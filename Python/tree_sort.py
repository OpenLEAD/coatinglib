from numpy import *
from scipy.spatial import KDTree
import time


approachrays = load('blade2.npz')
approachrays = approachrays['array']

T = KDTree(approachrays[:,0:3])

N = len(approachrays)
removeList = []
i=0
t=time.time()
while True:
    if i not in removeList:
        idx = T.query_ball_point(approachrays[i,0:3],r=0.01)
        idx = array(idx)
        idx = idx[idx>i]
        for j in idx:removeList.append(j)
        removeList.sort()
    if i%1000==0:
        print str(i)+'/'+str(N)
        print 'elapsed time= '+str(time.time()-t)
        t=time.time()
    i+=1
    if i==len(approachrays):break
approachrays = delete(approachrays,removeList,0)
