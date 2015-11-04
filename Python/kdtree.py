from numpy import *
from scipy.spatial import KDTree
import time


approachrays = load('blade2.npz')
approachrays = approachrays['array']
rays = approachrays[:,0:3]
T = KDTree(rays)

N = len(approachrays)
removeList = []
i=0
t=time.time()
while True:
    idx = T.query_ball_point(rays[i],r=0.01)
    idx = array(idx)
    idx = idx[idx>i]
    for j in idx:removeList.append(j)
    removeList = list(set(removeList))
    rays = delete(approachrays,removeList,0)
    rays = rays[:,0:3]
    if i%1000==0:
        print str(i)+'/'+str(N)+' and len(rays)='+str(len(rays))
        print 'elapsed time= '+str(time.time()-t)
        t=time.time()
    i+=1
    if i==len(rays):break
