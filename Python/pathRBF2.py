from numpy import *
import coating
from openravepy import *
from openravepy.misc import SpaceSamplerExtra
import time
import math
from time import gmtime, strftime
import RBF
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
dt = 3e-3 # parallel step
loops = 20
loopdistance = 1
#====================================================================================================================

# DIREITA E PARA BAIXO = 1
rR = RBF.rays
Tree = RBF.makeTree(rR)
Rn = 0

manipulabilityNUM = 0
globalManipPos = []
globalManipOri = []

def tangentOptm(ray,q0):
    tan = cross(ray[3:6],ray[0:3])
    tan = tan/sqrt(dot(tan,tan))
    res = coating.optmizeQ(robot,ikmodel,manip,ray,q0)

    angle_suc = isViable(res.x,ray[3:6])
    return tan, res.x, (res.success and angle_suc)

def solution(ray):
    _, iksolList, indexlist1 = coating.WorkspaceOnPose(pN, 0, [ray],robot,ikmodel,facevector,theta,coatingdistancetolerance)
    print 'solution: iksolList - ',array(iksolList).shape
    _, AlliksolList, _ = coating.AllExtraCoating2([ray],indexlist1,coatingdistance,numberofangles,tolerance,ikmodel,facevector,coatingdistancetolerance)
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
            _, iksolList, _ = coating.WorkspaceOnPose(pN, 0, [ray],robot,ikmodel,tempP,theta,coatingdistancetolerance)    
            if iksolList:
                return iksolList
    return []

def Qvector_backtrack(y,Q):
    for rev_y in reversed(y):
        res = coating.optmizeQ(robot,ikmodel,manip,rev_y,Q[-1])
        Q.append(res.x)
        if not res.success:
            print 'Qvector_backtrack error'
    return list(reversed(Q))

def Qvector(y,Q,sign):
    global handles
    suc=True
    while suc:
        tan, q, suc = tangentOptm(y[-1],Q[-1])
        if not coating.CheckDOFLimits(robot,q):
            #wait = input("deu bode, PRESS 1 ENTER TO CONTINUE.") 
            iksolList = fullSolution(y[-1])
            Q=[iksolList[0][0]]
            Q = Qvector_backtrack(y,Q)
            continue
        if suc:
            Q.append(q)
            p1 = curvepoint(y[-1][0:3]+sign*tan*dt)
            dv = array(RBF.df(p1))
            normdv = sqrt(dot(dv,dv))
            #print 'success:', res.success, ' fn:', polynomial_spline.fn4(xnew[0],xnew[1],xnew[2])/normdv
            n = dv/normdv
            P=concatenate((p1,n))   
            y.append(P)
            handles=plotPoint(P, handles,array((0,1,0)))
    return Q,y

def isViable(q0,norm):
    global manipulabilityNUM
    global globalManipPos
    global globalManipOri
    manipulability, manipjpos, manipjori = coating.manipulabilityS(manip)
    globalManipPos.append(manipjpos)
    globalManipOri.append(manipjori)
    #print 'manipjpos: ',manipjpos,', manipjori:', manipjori 

    manipjpos = 0.3*manipjpos + 0.7*manipulabilityNUM
    if manipulabilityNUM > manipjpos and manipulabilityNUM<=0.15:
        print "LOW MANIPULABILITY: ", manipjpos
        return False
    manipulabilityNUM = manipjpos
    
    if len(q0)>0:
        robot.SetDOFValues(q0,ikmodel.manip.GetArmIndices())
        T=manip.GetTransform()
        Rx = T[0:3,0]/sqrt(dot(T[0:3,0],T[0:3,0]))
        if dot(norm,Rx) >= math.cos(tolerance*math.pi/180):
            #print 'ANGLE: ', 180*math.acos(dot(norm,Rx))/math.pi
            return True
        else:
            print "ANGLE TOLERANCE FAIL: ", 180*math.acos(dot(norm,Rx))/math.pi
            return False
    else:
        return False

def drawParallel(ray,q0,sign):
    viable=isViable(q0,ray[3:6])
    Ql=[]; yL=[]; Qr=[]; yR=[]
    if viable:
        print 'drawParallel: is viable'
        QL=[q0]
        QR=[q0]
        y=[ray]
        #Arriscado - Pode caminhar em duas funcoes esquerda-direita
        if sign==1:
            Qr, yR = Qvector(y,QR,sign)
            #print 'sign 1: Qr -',array(Qr).shape,', yR -',array(yR).shape
            y=[ray]
            sign*=-1
            Ql, yL = Qvector(y,QL,sign)
            #print 'sign -1: Ql -',array(Ql).shape,', yL -',array(yL).shape
        else:
            Ql, yL = Qvector(y,QL,sign)
            y=[ray]
            sign*=-1
            Qr, yR = Qvector(y,QR,sign)
    else:
        sign*=-1
        print 'drawParallel: solution not found'
        while not viable:
            y=ray[0:3]
            tan, q0, viable = tangentOptm(ray,q0)
            y=curvepoint(y+sign*tan*dt)
            norm = RBF.df([y[0],y[1],y[2]])
            norm /= sqrt(dot(norm,norm))
            ray = array([y[0],y[1],y[2],norm[0],norm[1],norm[2]])
        print 'drawParallel: solution found'
        QL=[q0]
        QR=[q0]
        y=[ray]
        #Arriscado - Pode caminhar em duas funcoes esquerda-direita
        if sign==1:
            Qr, yR = Qvector(y,QR,sign)
        else:
            Ql, yL = Qvector(y,QL,sign)
        #return drawParallel(ray,q0,sign)
    
    print 'Ql -',array(list(reversed(Ql))).shape,', Qr -',array(Qr[1:]).shape
    
    if Ql and Qr[1:]:
        Q=concatenate((list(reversed(Ql)),Qr[1:]))
    elif Ql: Q=Ql
    else: Q=Qr
    return yR, yL, Q

def tangentd(ray,sign):
    P=ray
    x=P[0];y=P[1];z=P[2]
    dv = array(RBF.df([x,y,z]))
    n = dv/sqrt(dot(dv,dv))
    
    P=concatenate((P,n))

    tan = cross(n,[x,y,z])
    a = tan/sqrt(dot(tan,tan))
    tand = sign*cross(a,n)
    return tand

def tangent(ray,sign):
    x=ray[0];y=ray[1];z=ray[2]
    dv = array(RBF.df([x,y,z]))
    n = dv/sqrt(dot(dv,dv))
    tan = cross(n,[x,y,z])
    return tan/sqrt(dot(tan,tan))

def curvepoint(p0):
    tol=1e-4
    while True:
        df1 = RBF.df(p0); f1 = RBF.f(p0)
        df2 = array([2*p0[0],2*p0[1],2*p0[2]]); f2 = p0[0]**2+p0[1]**2+p0[2]**2-Rn**2
        df1df2 = dot(df1,df2); df1df1 = dot(df1,df1); df2df2 = dot(df2,df2)
        beta = (-f1+f2*df1df1/df1df2)/(df1df2-df1df1*df2df2/df1df2)
        alpha = (-f1+f2*df1df2/df2df2)/(df1df1-df1df2*df1df2/df2df2)
        dk = alpha*df1+beta*df2
        p1 = p0+dk
        dif = p1-p0
        if sqrt(dot(dif,dif))<tol:
            return p1
        else: p0=p1

def onecurvepoint(p0):
    tol=1e-4
    while True:
        df = RBF.df(p0)
        p1 = p0 - RBF.f(p0)*df/dot(df,df)
        dif = p1-p0
        if sqrt(dot(dif,dif))<tol:
            return p1
        else: p0=p1


def FindNextParallel(P0,sign):
    global Rn
    d = sqrt(P0[0]**2+P0[1]**2+P0[2]**2)
    R0=1.425
    n = math.ceil((d-R0)/0.003)
    Rn = R0+n*0.003
    y = curvepoint(P0)    
    norm = RBF.df([y[0],y[1],y[2]])
    norm /= sqrt(dot(norm,norm))
    y = array([y[0],y[1],y[2],norm[0],norm[1],norm[2]])        
    return y

def meridian2(P0,sign,q0):
    print 'meridian2'
    Q=[q0]
    P0=array([P0[0],P0[1],P0[2]])
    y = curvepoint(P0)
            
    norm = RBF.df(y)
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

def plotPoint(point, handles,color):
    handles.append(env.plot3(points=array(point)[0:3],pointsize=5,colors=color))
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
    Pd=FindNextParallel(P0,1) #1 para baixo
    q0=solution(Pd)
    return Pd, q0

def main():
    QALL = []
    global handles
    global Rn
    RBF.set_handles(handles,env)
    env.SetViewer('qtcoin')
    
    Pd, q0 = initialPoint()
    sign = 1
    q0=q0[0][0]
    robot.SetDOFValues(q0,ikmodel.manip.GetArmIndices())
    manipulabilityNUM = coating.manipulabilityDET(manip)
    
    for i in range(0,loops):
        yR, yL, Q = drawParallel(Pd,q0,sign)
        if sign==1:
            P0=yL[-1]
            qi=Q[0]
        else:
            P0=yR[-1]
            qi=Q[-1]
        QALL.extend(Q)
        Rn+=loopdistance*0.003
        Pd, Qd=meridian2(P0,1,qi)
        yM = getPointsfromQ(Qd)
        QALL.extend(Qd)
        handles=plotPoints(yM, handles,array((1,0,0)))
        
        sign*=-1
        q0=Qd[-1]
    return QALL

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

def getPointsfromQ(Q):
    points=[]
    for q in Q:
        robot.SetDOFValues(q,ikmodel.manip.GetArmIndices())
        T=manip.GetTransform()
        points.append(T[0:3,3])
    return array(points)  

if __name__ == '__main__':
    QALL = main()
    #savez_compressed('path/QALL/'+'t.npz', array=QALL)
    #wait = input("PRESS 1 ENTER TO CONTINUE.")
    #coating.robotPath2(QALL, 0.005,robot,ikmodel)
#coating.robotPath2(Q, 0.005,robot,ikmodel)
#env.SetViewer('qtcoin')
#handles=[]
#handles.append(env.plot3(points=array(yL)[:,0:3],pointsize=5,colors=array((0,0,0))))
#handles.append(env.plot3(points=yR[:,0:3],pointsize=5,colors=array((0,0,0))))
