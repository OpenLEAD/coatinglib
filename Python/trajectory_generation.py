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
target = env.GetBodies()[0]
handles=[]

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

dt = 3e-3 # parallel step
loops = 1
loopdistance = 1
#====================================================================================================================

# DIREITA E PARA BAIXO = 1
rR = RBF.rays
Tree = RBF.makeTree(rR)
Rn = 0

def tangent(ray):
    tan = cross(ray[3:6],ray[0:3])
    tan = tan/sqrt(dot(tan,tan))   
    return tan

def drawParallel(ray):
    global handles
    y = [ray]
    while True:
        tan = tangent(y[-1])
        #print 'tan= ', tan[0],tan[1],tan[2]
        p1 = curvepoint(y[-1][0:3]-tan*dt)
        dv = array(RBF.df(p1))
        normdv = sqrt(dot(dv,dv))
        n = dv/normdv
        P=concatenate((p1,n))   
        y.append(P)
        handles=plotPoint(P, handles,array((0,1,0)))
    return y

def curvepoint(p0):
    tol=1e-4
    while True: 
        df1 = RBF.df(p0); f1 = RBF.f(p0)
        df2 = array([2*p0[0],2*p0[1],2*p0[2]]); f2 = p0[0]**2+p0[1]**2+p0[2]**2-Rn**2
        df1df2 = dot(df1,df2); df1df1 = dot(df1,df1); df2df2 = dot(df2,df2)
        beta = (-f1+f2*df1df1/df1df2)/(df1df2-df1df1*df2df2/df1df2)
        alpha = (-f1+f2*df1df2/df2df2)/(df1df1-df1df2*df1df2/df2df2)
        dk = alpha*df1+beta*df2
        #print 'preso= ', sqrt(dot(dk,dk))
        p1 = p0+dk
        if sqrt(dot(dk,dk))<tol:
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
    #soma=[0,0,0]
    #for i in rRp:soma=soma+i
    #P0 = soma/len(rRp)
    P0 = rRp[argmin(rRp[:,1])]
    P0 = nearestSample(P0,rRp)
    Pd=FindNextParallel(P0,1) #1 para baixo
    return Pd

def main():
    global handles
    global Rn
    RBF.set_handles(handles,env)
    env.SetViewer('qtcoin')
    
    Pd = initialPoint()
    for i in range(0,loops):
        y = drawParallel(Pd)
        P0=y[-1]
        Rn+=loopdistance*0.003
        #Pd, Qd=meridian2(P0,1,qi)
    return y

if __name__ == '__main__':
    y = main()
    #savez_compressed('path/QALL/'+'t.npz', array=QALL)
    #wait = input("PRESS 1 ENTER TO CONTINUE.")
    #coating.robotPath2(QALL, 0.005,robot,ikmodel)
#coating.robotPath2(Q, 0.005,robot,ikmodel)
#env.SetViewer('qtcoin')
#handles=[]
#handles.append(env.plot3(points=array(yL)[:,0:3],pointsize=5,colors=array((0,0,0))))
#handles.append(env.plot3(points=yR[:,0:3],pointsize=5,colors=array((0,0,0))))
