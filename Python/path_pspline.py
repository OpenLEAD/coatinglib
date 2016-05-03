from numpy import *
import min_distance
import coating
from openravepy import *
from openravepy.misc import SpaceSamplerExtra
import time
import math
from time import gmtime, strftime
import polynomial_spline
#====================================================================================================================
env=Environment()
env.Load("/home/renan/workspace/coatinglib/Turbina/env_mh12_0_16.xml")
robot = env.GetRobots()[0]
target = env.GetBodies()[0]
manip = robot.GetActiveManipulator()
handles=[]

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
T = array([[1,0,0,0],[0,cos(alpha),-sin(alpha),0],
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
pN = concatenate((pN,normal))
QList = []
#====================================================================================================================

# DIREITA E PARA BAIXO = 1
rR = polynomial_spline.rays
Tree = polynomial_spline.makeTree(rR)

def tangentOptm(ray,q0):
    tan = cross(ray[3:6],ray[0:3])
    tan *= (40.0/60)/sqrt(dot(tan,tan))
    angle_suc = True
    limits = True
    #print 'ray - ', ray[0:3]
    res = coating.optmizeQ(robot,ikmodel,manip,ray,q0)
    manipulability = coating.manipulabilityDET(manip)>=1e-3
    if res.success and not isViable(res.x,ray[3:6]):
        angle_suc=False    
    return tan, res.x, (res.success and angle_suc and manipulability)

def solution(ray):
    reachableRays, iksolList, indexlist1 = coating.WorkspaceOnPose(pN, 0, [ray],robot,ikmodel,facevector,theta,coatingdistancetolerance)
    print 'solution: iksolList - ',array(iksolList).shape
    AllreachableRays, AlliksolList, indexlist2 = coating.AllExtraCoating2([ray],indexlist1,coatingdistance,numberofangles,tolerance,ikmodel,facevector,coatingdistancetolerance)
    print 'solution: AlliksolList - ',array(AlliksolList).shape
    if iksolList: return iksolList
    elif AlliksolList: return AlliksolList
    else: return False

def fullSolution(ray):
    angles = arange(0,tolerance,0.1)
    numberofangles = 10
    for angle in angles:
        angle=1.0*pi*angle/180
        Rv3tol = coating.RotationCalc2([0,0,1],angle)
        p = dot(facevector,transpose(Rv3tol))
        k = 1.0*2*pi/numberofangles
        for i in range(0,numberofangles):
            alfa = k*i
            Rv1alfa = coating.RotationCalc2(facevector,alfa)
            tempP = dot(p,transpose(Rv1alfa))
            reachableRays, iksolList, indexlist1 = coating.WorkspaceOnPose(pN, 0, [ray],robot,ikmodel,tempP,theta,coatingdistancetolerance)    
            if iksolList:
                return iksolList
    return []

def Qvector(y,Q,dt,sign):
    suc=True
    while suc:
        tan, q, suc = tangentOptm(y[-1],Q[-1])
        if polynomial_spline.fn4(y[-1][0],y[-1][1],y[-1][2]) > 1e-6 :
            print polynomial_spline.fn4(y[-1][0],y[-1][1],y[-1][2])
        if not coating.CheckDOFLimits(robot,q):
                return [], y[-1]
        if suc:
            Q.append(q)
            xold = y[-1][0:3]
            pnew = xold + sign*tan*dt
            tol=1e-5
            xnew = xold
            while dot(xnew-xold,xnew-xold)<=tol: # ==0
                res = coating.optmizeTan(xold, pnew, tol)
                tol*=0.1
                xnew = res.x
            dv = array(polynomial_spline.dfn4(xnew[0],xnew[1],xnew[2]))
            n = dv/sqrt(dot(dv,dv))
            P=concatenate((xnew,n))   
            y.append(P)
            
    return Q,y

def isViable(q0,norm):
    robot.SetDOFValues(q0,ikmodel.manip.GetArmIndices())
    T=manip.GetTransform()
    Rx = T[0:3,0]/sqrt(dot(T[0:3,0],T[0:3,0]))
    return dot(norm,Rx) >= math.cos(tolerance*math.pi/180)

def drawParallel2(ray,q0,sign):
    viable=isViable(q0,ray[3:6])
    dt = 1e-3
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
        if not Qr:
            iksolList = fullSolution(yR)
            return drawParallel2(yR,iksolList[0][0],sign*-1)
        if not Ql:
            iksolList = fullSolution(yL)
            return drawParallel2(yL,iksolList[0][0],sign*-1)
    else:
        sign*=-1
        print 'drawParallel2: solution not found'
        tol = 1e-6
        while not viable:
            y=ray[0:3]
            tan, q0, viable = tangentOptm(ray,q0)
            tan *= sign*dt
            pnew = y+tan
            res = coating.optmizeTan(y, pnew, tol)
            if dot(res.x-y,res.x-y)==0:
                tol*=0.1
            y=res.x

            norm = polynomial_spline.dfn4(y[0],y[1],y[2])
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
    dt = 1e-3
    if sol:
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

            
            norm = polynomial_spline.dfn4(y[0],y[1],y[2])
            norm /= sqrt(dot(norm,norm))       
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
    dv = array(polynomial_spline.dfn4(x,y,z))
    n = dv/sqrt(dot(dv,dv))
    
    P=concatenate((P,n))

    tan = cross(n,[x,y,z])
    a = tan/sqrt(dot(tan,tan))
    tand = sign*cross(a,n)
    return tand

def tangent(ray,sign):
    P=ray
    x=P[0];y=P[1];z=P[2]
    dv = array(polynomial_spline.dfn4(x,y,z))
    n = dv/sqrt(dot(dv,dv))
    
    P=concatenate((P,n))

    tan = cross(n,[x,y,z])
    a = tan/sqrt(dot(tan,tan))
    return a

def meridian(P0,Rn,sign):
    dt = 1e-4
    y=array([float(P0[0]),float(P0[1]),float(P0[2])])
    d = sqrt(P0[0]**2+P0[1]**2+P0[2]**2)
    R0=1.425
    n = math.ceil((d-R0)/0.003)
    Rn = R0+n*0.003
    dif = d-Rn
    S=(dif>0)
    tol=1e-6
    while abs(dif)>1e-4:
        tand = tangentd(y,sign)*dt
        pnew = y+tand
        res = coating.optmizeTan(y, pnew, tol)
        
        if dot(res.x-y,res.x-y)==0:tol*=0.1
        y = res.x
        d = sqrt(res.x[0]**2+res.x[1]**2+res.x[2]**2)
        dif = d-Rn
        if S!=(dif>0):
            sign*=-1
            dt*=0.5
            S=(dif>0)
    
    norm = polynomial_spline.dfn4(y[0],y[1],y[2])
    norm /= sqrt(dot(norm,norm))
    y = array([y[0],y[1],y[2],norm[0],norm[1],norm[2]])        
    return y

def meridian2(P0,Rn,sign,q0):
    print 'meridian2'
    Q=[q0]
    dt = 1e-4
    y=array([float(P0[0]),float(P0[1]),float(P0[2])])
    d = sqrt(P0[0]**2+P0[1]**2+P0[2]**2)
    dif = d-Rn
    S=(dif>0)
    tol=1e-6
    notstop = True
    while abs(dif)>1e-4:
        tand = tangentd(y,sign)*dt
        pnew = y+tand
        res = coating.optmizeTan(y, pnew,tol)
        if dot(res.x-y,res.x-y)==0:tol*=0.1
        y=res.x

        if notstop:
            norm = polynomial_spline.dfn4(y[0],y[1],y[2])
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

            
    norm = polynomial_spline.dfn4(y[0],y[1],y[2])
    norm /= sqrt(dot(norm,norm))
    ray = array([y[0],y[1],y[2],norm[0],norm[1],norm[2]])  

    resQ = coating.optmizeQ(robot,ikmodel,manip,ray,Q[-1])
    if resQ.success:
        Q.append(resQ.x)
    else:
        print 'Meridian Error'
        
    return ray, Q

def nearestSample(P,points):
    _,idx = Tree.query(P)
    proximo=points[idx,0:3]
    return proximo
   
def plotPoints(points, handles,color):
    handles.append(env.plot3(points=array(points)[:,0:3],pointsize=5,colors=color))
    return handles

def nextLine(P0):
    d0 = sqrt(P0[0]**2+P0[1]**2+P0[2]**2)
    R0=1.425
    n = math.ceil((d0-R0)/0.003)
    Rn = R0+n*0.003
    return Rn

def initialPoint():
    rRp = rR[:,0:3]
    soma=[0,0,0]
    for i in rRp:soma=soma+i
    P0 = soma/len(rRp)
    P0 = nearestSample(P0,rRp)
    
    Pd=meridian(P0,Rn,1)
    return Pd

def main():
    env.SetViewer('qtcoin')
    Pd = initialPoint()
    yR, yL, Q = drawParallel(Pd,1)
    t=strftime("%d%b%Y_%H_%M_%S", gmtime())
    savez_compressed('path/yR/'+'_'+t+'.npz', array=yR)
    savez_compressed('path/yL/'+'_'+t+'.npz', array=yL)
    savez_compressed('path/Q/'+'_'+t+'.npz', array=Q)
    
    return yR, yL, Q

def getPointsfromQ(Q):
    points=[]
    for q in Q:
        robot.SetDOFValues(q,ikmodel.manip.GetArmIndices())
        T=manip.GetTransform()
        points.append(T[0:3,3])
    return array(points)    


def main2():
    global handles
    polynomial_spline.set_handles(handles,env)
    #env.SetViewer('qtcoin')
    yR, yL, Q = main()
    handles=plotPoints(yL, handles,array((0,0,0)))
    handles=plotPoints(yR, handles,array((0,0,0)))
    P0=yL[-1]
    Rn=nextLine(P0)
    Rn+=33*0.003
    Pd, Qd=meridian2(P0,Rn,1,Q[0])
    yM = getPointsfromQ(Qd)
    yR2, yL2, Q2 = drawParallel2(Pd,Qd[-1],-1)
    handles=plotPoints(yM, handles,array((0,0,0)))
    handles=plotPoints(yL2, handles,array((0,0,0)))
    handles=plotPoints(yR2, handles,array((0,0,0)))

    P0=yR2[-1]
    Rn=nextLine(P0)
    Rn+=33*0.003
    Pd, Qd=meridian2(P0,Rn,1,Q[-1])
    yM = getPointsfromQ(Qd)
    yR3, yL3, Q3 = drawParallel2(Pd,Qd[-1],1)
    handles=plotPoints(yM, handles,array((0,0,0)))
    handles=plotPoints(yL3, handles,array((0,0,0)))
    handles=plotPoints(yR3, handles,array((0,0,0)))
    
    return yR3, yL3, Q3,yR2, yL2, Q2, yR, yL, Q

def realCoatedPoints(Q):
    global handles
    newY=[]
    for q in Q:
        robot.SetDOFValues(q,ikmodel.manip.GetArmIndices())
        T=manip.GetTransform()
        Rx = array(T[0:3,0]/sqrt(dot(T[0:3,0],T[0:3,0])))
        p = T[0:3,3]
        newY.append(p-0.23*Rx)
    handles=plotPoints(newY, handles,array((0,0,1)))
    return handles

#coating.robotPath2(Q, 0.005,robot,ikmodel)
#env.SetViewer('qtcoin')
#handles=[]
#handles.append(env.plot3(points=array(yL)[:,0:3],pointsize=5,colors=array((0,0,0))))
#handles.append(env.plot3(points=yR[:,0:3],pointsize=5,colors=array((0,0,0))))
