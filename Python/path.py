from numpy import *
from moving_LS import *
import matplotlib.pyplot as plt
import min_distance
import coating
from openravepy import *
from openravepy.misc import SpaceSamplerExtra
import time
import math
import moving_LS

#====================================================================================================================
env=Environment()
env.Load("/home/renan/workspace/coatinglib/Turbina/env_mh12_0_16.xml")
robot = env.GetRobots()[0]
target = env.GetBodies()[0]
manip = robot.GetActiveManipulator()

ikmodel = databases.inversekinematics.InverseKinematicsModel(robot=robot,iktype=IkParameterization.Type.Transform6D)
if not ikmodel.load():
    ikmodel.autogenerate()

# PARAMETERS
facevector = [1,0,0]
theta = [0,0,0]
coatingdistance = 0.23 # coating distance
coatingdistancetolerance = 0.01
numberofangles = 8 # degree step
tolerance = 30 # degrees
alpha = 1.0*pi/180; #degree blade step
BladePosition = 0

# Initial T
Ti = []
for body in env.GetBodies():
    Ti.append(body.GetTransform())

alpha = 1.0*BladePosition*pi/180
T = numpy.array([[1,0,0,0],[0,cos(alpha),-sin(alpha),0],
                 [0,sin(alpha),cos(alpha),0],[0,0,0,1]])

for ibody in range(0,len(env.GetBodies())):
    if ibody!=5:
        env.GetBodies()[ibody].SetTransform(dot(T,Ti[ibody]))

gbase_position=array([[ -2.14387087e+00,  -3.22000000e+00,  -9.97857391e-01],
       [ -1.98593217e+00,  -3.22000000e+00,  -6.64834872e-01],
       [ -1.82799346e+00,  -3.22000000e+00,  -3.31812352e-01],
       [ -1.67005475e+00,  -3.22000000e+00,   1.21016736e-03],
       [ -1.51211605e+00,  -3.22000000e+00,   3.34232687e-01],
       [ -1.35417734e+00,  -3.22000000e+00,   6.67255206e-01],
       [ -1.19623863e+00,  -3.22000000e+00,   1.00027773e+00]])

pN=gbase_position[4]
normal = [-1,0,0]
pN = numpy.concatenate((pN,normal))
QList = []
#====================================================================================================================

# DIREITA E PARA BAIXO = 1
        
reachableRays = load('coated_points/pos4_black.npz')
reachableRays  = reachableRays['array']
AllreachableRays = load('coated_points/pos4_blue.npz')
AllreachableRays  = AllreachableRays['array']
rR = concatenate((reachableRays,AllreachableRays))

v = array([  6.11242244e-02,   2.89653535e-01,   7.94314249e-01,
        -3.20708568e-02,   5.11681485e-02,   5.32913320e-01,
        -1.78331728e-01,   4.31527747e-01,  -1.08549881e-02,
         3.21256789e-01,  -7.03089123e-02,   7.22629763e-02,
         7.97448256e-02,  -8.52081772e-04,   7.66169966e-03,
        -1.13579708e-01,  -6.03034682e-01,  -6.98038219e-01,
        -1.89876807e-01,   1.56603098e-01,  -5.04314372e-01,
        -3.11011260e-02,  -7.67252299e-02,  -1.48860963e-01,
        -3.31170392e-02,   5.83401508e-01,   5.32945098e-01,
         3.44663122e-01,   9.20029194e-01,  -3.58423952e-01,
         2.63214686e-01,  -9.43235599e-02,  -9.33102549e-02,
         3.41146158e-01,  -1.94242634e-02])

def vector4(x,y,z):
    return [1, z, z**2, z**3, z**4, y, y*z, y*z**2, y*z**3, y**2, y**2*z, y**2*z**2, y**3, y**3*z, y**4, x, x*z, x*z**2, x*z**3, x*y, x*y*z, x*y*z**2, x*y**2, x*y**2*z, x*y**3, x**2, x**2*z, x**2*z**2, x**2*y, x**2*y*z, x**2*y**2, x**3, x**3*z, x**3*y, x**4]

def dvector4(x,y,z):
    dx = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, z, z**2, z**3, y, y*z, y*z**2, y**2, y**2*z, y**3, 2*x, 2*x*z, 2*x*z**2, 2*x*y, 2*x*y*z, 2*x*y**2, 3*x**2, 3*x**2*z, 3*x**2*y, 4*x**3]
    dy = [0, 0, 0, 0, 0, 1, z, z**2, z**3, 2*y, 2*y*z, 2*y*z**2, 3*y**2, 3*y**2*z, 4*y**3, 0, 0, 0, 0, x, x*z, x*z**2, 2*x*y, 2*x*y*z, 3*x*y**2, 0, 0, 0, x**2, x**2*z, 2*x**2*y, 0, 0, x**3, 0]
    dz = [0, 1, 2*z, 3*z**2, 4*z**3, 0, y, 2*y*z, 3*y*z**2, 0, y**2, 2*y**2*z, 0, y**3, 0, 0, x, 2*x*z, 3*x*z**2, 0, x*y, 2*x*y*z, 0, x*y**2, 0, 0, x**2, 2*x**2*z, 0, x**2*y, 0, 0, x**3, 0, 0]
    a=[dx,dy,dz]
    return a

def fn4(x,y,z):
    v, _, _, _, _ = polynomial_surface(array([x,y,z]),rR)
    return dot(v,vector4(x,y,z))
#    return dot(vector4(x,y,z),transpose(v))

def dfn4(x,y,z):
    return dpolynomial(array([x,y,z]),rR)
    #dv = dvector4(x,y,z)
    #return [dot(dv[0],transpose(v)),dot(dv[1],transpose(v)),dot(dv[2],transpose(v))]

def getmean(points):
    m=[0,0,0,0,0,0]
    for p in points:m=m+p
    m=m/len(points)

    pmean=[0,0,0,0,0,0]
    pmin=1000
    for p in points:
        d=dot(m[0:3]-p[0:3],m[0:3]-p[0:3])
        if pmin>d:
            pmin=d
            pmean=p
    return pmean

def tangentOptm(ray,q0):
    tan = cross(ray[3:6],ray[0:3])
    tan *= (40.0/60)/sqrt(dot(tan,tan))
 
    angle_suc = True
    res = coating.optmizeQ(robot,ikmodel,manip,ray,q0)

##    robot.SetDOFValues(res.x,ikmodel.manip.GetArmIndices())
##    T=manip.GetTransform()
##    Rx = T[0:3,0]/sqrt(dot(T[0:3,0],T[0:3,0]))
   # print 'tangentOptm: res.success -',res.success , ' angle -', 180*math.acos(dot(ray[3:6],Rx))/math.pi
    if res.success and not isViable(res.x,ray[3:6]):
        angle_suc=False
       

    return tan, res.x, (res.success and angle_suc)

def HasSol(ray):
    reachableRays, iksolList, indexlist1 = coating.WorkspaceOnPose(pN, 0, [ray],robot,ikmodel,facevector,theta,coatingdistancetolerance)
    AllreachableRayschableRays, AlliksolList, indexlist2 = coating.AllExtraCoating2([ray],indexlist1,coatingdistance,numberofangles,tolerance,ikmodel,facevector,coatingdistancetolerance)
    if AlliksolList: return True
    else: return False

def solution(ray):
    print 'solution: tolerance - ', tolerance
    reachableRays, iksolList, indexlist1 = coating.WorkspaceOnPose(pN, 0, [ray],robot,ikmodel,facevector,theta,coatingdistancetolerance)
    print 'solution: iksolList - ',array(iksolList).shape
    AllreachableRays, AlliksolList, indexlist2 = coating.AllExtraCoating2([ray],indexlist1,coatingdistance,numberofangles,tolerance,ikmodel,facevector,coatingdistancetolerance)
    print 'solution: AlliksolList - ',array(AlliksolList).shape
    print AlliksolList
    if iksolList: return iksolList
    elif AlliksolList: return AlliksolList
    else: return False
    
def Qvector(y,Q,dt,sign):
    suc=True
    while suc:
        tan, q, suc = tangentOptm(y[-1],Q[-1])
        #print 'Qvec: tan -',array(tan).shape,', Q -',array(Q).shape,', suc -',suc

        if suc:
            Q.append(q)
            xold = y[-1][0:3]
            pnew = xold + sign*tan*dt
            tol=1e-5
            xnew = xold
            while dot(xnew-xold,xnew-xold)==0:
                res = coating.optmizeTan(xold, pnew, tol)
                tol*=0.1
                xnew = res.x

            dv = array(dfn4(xnew[0],xnew[1],xnew[2]))
            n = dv/sqrt(dot(dv,dv))
            
            P=concatenate((xnew,n))
            
            y.append(P)
    return Q,y

def isViable(q0,norm):
    robot.SetDOFValues(q0,ikmodel.manip.GetArmIndices())
    T=manip.GetTransform()
    Rx = T[0:3,0]/sqrt(dot(T[0:3,0],T[0:3,0]))
    return dot(norm,Rx) >= math.cos(30*math.pi/180)

def drawParallel2(ray,q0,sign):
    print 'drawparallel2'
    viable=isViable(q0,ray[3:6])
    dt = 0.0015
    if viable:
        QL=[q0]
        QR=[q0]
        y=[ray]
        if sign==1:
            Qr, yR = Qvector(y,QR,dt,sign)
            #print 'sign 1: Qr -',array(Qr).shape,', yR -',array(yR).shape
            y=[ray]
            sign*=-1
            Ql, yL = Qvector(y,QL,dt,sign)
            #print 'sign -1: Ql -',array(Ql).shape,', yL -',array(yL).shape

        else:
            Ql, yL = Qvector(y,QL,dt,sign)
            y=[ray]
            sign*=-1
            Qr, yR = Qvector(y,QR,dt,sign)
    else:
        sign*=-1
        print 'drawParallel2: solution not found'
        tol = 1e-6
        while not viable:
            y=ray[0:3]
            tan, q0, viable = tangentOptm(ray,q0)
            tan *= sign*dt
            pnew = y+tan
            res = coating.optmizeTan(y, pnew, v,tol)
            if dot(res.x-y,res.x-y)==0:
                tol*=0.1
            y=res.x

            norm = dfn4(y[0],y[1],y[2])
            norm /= sqrt(dot(norm,norm))
            ray = array([y[0],y[1],y[2],norm[0],norm[1],norm[2]])
        return drawParallel2(ray,q0,sign)
    
    print 'Ql -',array(list(reversed(Ql))).shape,', Qr -',array(Qr[1:]).shape
    
    if Ql and Qr[1:]:
        Q=concatenate((list(reversed(Ql)),Qr[1:]))
    elif Ql: Q=Ql
    else: Q=Qr
    return yR, yL, Q

def drawParallel(ray,sign):
    sol=solution(ray)
    dt = 0.0015
    if sol:
        print len(sol[0])
        q0=sol[0][0]
        QL=[q0]
        QR=[q0]
        y=[ray]
        if sign==1:
            Qr, yR = Qvector(y,QR,dt,sign)
            #print 'sign 1: Qr -',array(Qr).shape,', yR -',array(yR).shape
            y=[ray]
            sign*=-1
            Ql, yL = Qvector(y,QL,dt,sign)
            #print 'sign -1: Ql -',array(Ql).shape,', yL -',array(yL).shape

        else:
            Ql, yL = Qvector(y,QL,dt,sign)
            y=[ray]
            sign*=-1
            Qr, yR = Qvector(y,QR,dt,sign)
    else:
        sign*=-1
        print 'drawParallel: solution not found'
        tol = 1e-6
        while not solution(ray):
            y=ray[0:3]
            tan = tangent(y,sign)*dt
            pnew = y+tan
            res = coating.optmizeTan(y, pnew, v,tol)
            if dot(res.x-y,res.x-y)==0:tol*=0.1
            y=res.x

            
            norm = dfn4(y[0],y[1],y[2])
            norm /= sqrt(dot(norm,norm))       

##            a=dvector4(y[0],y[1],y[2])
##            norm=[dot(a[0],transpose(v)),dot(a[1],transpose(v)),dot(a[2],transpose(v))]
            ray = [y[0],y[1],y[2],norm[0],norm[1],norm[2]]
        return drawParallel(ray,sign)
    
    print 'Ql -',array(list(reversed(Ql))).shape,', Qr -',array(Qr[1:]).shape
    
    if Ql and Qr[1:]:
        Q=concatenate((list(reversed(Ql)),Qr[1:]))
    elif Ql: Q=Ql
    else: Q=Qr
    return yR, yL, Q

def tangentd(ray,sign):
    P=ray
    x=P[0];y=P[1];z=P[2]
    dv = array(dfn4(x,y,z))
    n = dv/sqrt(dot(dv,dv))
    
    P=concatenate((P,n))

    tan = cross(n,[x,y,z])
    a = tan/sqrt(dot(tan,tan))
    tand = sign*cross(a,n)
    return tand

def tangent(ray,sign):
    P=ray
    x=P[0];y=P[1];z=P[2]
    dv = array(dfn4(x,y,z))
    n = dv/sqrt(dot(dv,dv))
    
    P=concatenate((P,n))

    tan = cross(n,[x,y,z])
    a = tan/sqrt(dot(tan,tan))
    return a

def meridian(P0,Rn,sign):
    dt = 0.001
    y=array([float(P0[0]),float(P0[1]),float(P0[2])])
    d = sqrt(P0[0]**2+P0[1]**2+P0[2]**2)
    dif = d-Rn
    S=(dif>0)
    tol=1e-6
    while abs(dif)>1e-4:
        tand = tangentd(y,sign)*dt
        pnew = y+tand
        res = coating.optmizeTan(y, pnew, tol)
        if dot(res.x-y,res.x-y)==0:tol*=0.1
        y=res.x
        d = sqrt(res.x[0]**2+res.x[1]**2+res.x[2]**2)
        dif = d-Rn
        if S!=(dif>0):
            sign*=-1
            dt*=0.5
            S=(dif>0)

    norm = dfn4(y[0],y[1],y[2])
    norm /= sqrt(dot(norm,norm))
    y = array([y[0],y[1],y[2],norm[0],norm[1],norm[2]])        
    return y

def meridian2(P0,Rn,sign,q0):
    print 'meridian2'
    Q=[q0]
    dt = 0.001
    y=array([float(P0[0]),float(P0[1]),float(P0[2])])
    d = sqrt(P0[0]**2+P0[1]**2+P0[2]**2)
    dif = d-Rn
    S=(dif>0)
    tol=1e-6
    notstop = True
    while abs(dif)>1e-4:
        tand = tangentd(y,sign)*dt
        pnew = y+tand
        res = coating.optmizeTan(y, pnew, v,tol)
        if dot(res.x-y,res.x-y)==0:
            tol*=0.1
        y=res.x

        if notstop:
            norm = dfn4(y[0],y[1],y[2])
            norm /= sqrt(dot(norm,norm))
            ray = array([y[0],y[1],y[2],norm[0],norm[1],norm[2]])

            resQ = coating.optmizeQ(robot,ikmodel,manip,ray,Q[-1])
            if resQ.success:
                Q.append(resQ.x)
            else:
                print 'Meridian Error'
        
        d = sqrt(res.x[0]**2+res.x[1]**2+res.x[2]**2)
        dif = d-Rn
        if S!=(dif>0):
            notstop = False
            sign*=-1
            dt*=0.5
            S=(dif>0)

            
    norm = dfn4(y[0],y[1],y[2])
    norm /= sqrt(dot(norm,norm))
    ray = array([y[0],y[1],y[2],norm[0],norm[1],norm[2]])  

    resQ = coating.optmizeQ(robot,ikmodel,manip,ray,Q[-1])
    if resQ.success:
        Q.append(resQ.x)
    else:
        print 'Meridian Error'
        
    return ray, Q

def proximityInCloud(P,points):
    dmin=1000
    proximo=[0,0,0]
    for point in points:
        dP = P-point
        d = dot(dP,dP)
        if d<dmin:
            dmin=d
            proximo=point
    return proximo        

def nearInSurface(P,points):
    proximo = proximityInCloud(P,points)
    P = min_distance.minPoint(proximo)
    return P
   
def plotPoints(points, handles,color):
    handles.append(env.plot3(points=array(points)[:,0:3],pointsize=5,colors=color))
    return handles

def nextLine(P0):
    d0 = sqrt(P0[0]**2+P0[1]**2+P0[2]**2)
    R0=1.425
    n = math.ceil((d0-R0)/0.003)
    Rn = R0+n*0.003
    return Rn

def main():
    rRp = rR[:,0:3]
    soma=[0,0,0]
    for i in rRp:soma=soma+i
    P0 = soma/len(rRp)
    P0 = nearInSurface(P0,rRp)
    d0 = sqrt(P0[0]**2+P0[1]**2+P0[2]**2)
    R0=1.425
    n = math.ceil((d0-R0)/0.003)
    Rn = R0+n*0.003
    
    Pd=meridian(P0,Rn,1)

    yR, yL, Q = drawParallel(Pd,1)

    
    return yR, yL, Q

def getPointsfromQ(Q):
    points=[]
    for q in Q:
        robot.SetDOFValues(q,ikmodel.manip.GetArmIndices())
        T=manip.GetTransform()
        points.append(T[0:3,3])
    return array(points)    

handles=[]
def main2():
    global handles
    env.SetViewer('qtcoin')
    yR, yL, Q = main()
    handles=plotPoints(yL, handles,array((0,0,0)))
    handles=plotPoints(yR, handles,array((0,0,0)))
    P0=yR[-1]
    Rn=nextLine(P0)
    Rn+=33*0.003
    Pd, Qd=meridian2(P0,Rn,1,Q[-1])
    yM = getPointsfromQ(Qd)
    yR2, yL2, Q2 = drawParallel2(Pd,Qd[-1],1)
    handles=plotPoints(yM, handles,array((0,1,0)))
    handles=plotPoints(yL2, handles,array((1,0,0)))
    handles=plotPoints(yR2, handles,array((1,0,0)))

##    P0=yR2[-1]
##    Rn=nextLine(P0)
##    Rn+=33*0.003
##    Pd, Qd=meridian2(P0,Rn,1,Q[0])
##    yM = getPointsfromQ(Qd)
##    yR2, yL2, Q2 = drawParallel2(Pd,Qd[-1],-1)
##    handles=plotPoints(yM, handles,array((0,1,0)))
##    handles=plotPoints(yL2, handles,array((0,0,1)))
##    handles=plotPoints(yR2, handles,array((0,0,1)))
    
    return yR2, yL2, Q2, yR, yL, Q
#coating.robotPath2(Q, 0.005,robot,ikmodel)
#env.SetViewer('qtcoin')
#handles=[]
#handles.append(env.plot3(points=array(yL)[:,0:3],pointsize=5,colors=array((0,0,0))))
#handles.append(env.plot3(points=yR[:,0:3],pointsize=5,colors=array((0,0,0))))
