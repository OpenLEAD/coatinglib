from numpy import *
from scipy.linalg import solve

trajectories = load('trajectory/Y.npz')
trajectories = trajectories['array']

for trajectory in trajectories:
    for index in range(0,len(trajectory)):
        p0 = trajectory[index-1][0:3]
        p1 = trajectory[index][0:3]
        p2 = trajectory[index+1][0:3]
        K = zeros((3,3))
        d = zeros(3)
        for i in [p0,p1,p2]:
            K[0,0] += i[0]**2
            K[0,1] += i[0]*i[1]
            K[0,2] += i[0]
            K[1,1] += i[1]**2
            K[1,2] += i[1]
            d[0] += -i[0]*i[2]
            d[1] += -i[2]*i[1]
            d[2] += -i[2]
        K[1,0]=K[0,1]    
        K[2,0]=K[0,2]
        K[2,1]=K[1,2]
        K[2,2]=3
        a=solve(K,d)
        c = a[2]
        a=array([a[0],a[1],1])
        projectedPoints = []
        for i in [p0,p1,p2]:
            dist = dot(i,a)+c
            projectedPoints.append(i-dist*a)
        projectedPoints=array(projectedPoints)
        t = projectedPoints[1]-projectedPoints[0]
        u = projectedPoints[2]-projectedPoints[0]
        v = projectedPoints[2]-projectedPoints[1]
        w = cross(t,u)
        iwsl2 = 1.0/(2*sqrt(dot(w,w)))
        radius = sqrt(dot(t,t)*dot(u,u)*dot(v,v)*iwsl2)
            
