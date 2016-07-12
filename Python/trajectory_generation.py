#!/usr/bin/env python
# Given the blade model (RBF) and the samples,
# this script generate the robot trajectories. The trajectories are the
# intersection between the spheres (radius 1.425 to 3.770) with the RBF model of
# the blade. The algorithm follows the marching method, available in
# http://www.mathematik.tu-darmstadt.de/~ehartmann/cdgen0104.pdf page 94
# intersection surface - surface.

# TODO: parallel computation of the trajectories as they are independent
# of the radius of the sphere
import coating
from numpy import *
from openravepy import *
from openravepy.misc import SpaceSamplerExtra
import math
import RBF
#====================================================================================================================
env=Environment()
env.Load("/home/renan/workspace/coatinglib/Turbina/env_mh12_0_16.xml")
target = env.GetBodies()[0]
handles=[]

dt = 3e-3 # parallel step
loops = 1
loopdistance = 1
rR = RBF.rays
Tree = RBF.makeTree(rR)
Rn = 0
#====================================================================================================================

def tangent(ray):
    tan = cross(ray[3:6],ray[0:3])
    tan = tan/sqrt(dot(tan,tan))   
    return tan

def drawParallel(Y,Pd):
    global handles
    theta0 = 180*math.atan2(-Pd[2],Pd[0])/math.pi
    counter = 0
    y = [Pd]       
    while True:
        tan = tangent(y[-1])
        p1 = curvepoint(y[-1][0:3]-tan*dt)
        dv = array(RBF.df(p1))
        normdv = sqrt(dot(dv,dv))
        n = dv/normdv
        P=concatenate((p1,n))

        r = math.sqrt(P[0]**2+P[1]**2+P[2]**2)
        theta = 180*math.atan2(-P[2],P[0])/math.pi
        if counter==0:
            if theta>theta0:
                counter+=1
        else:
            if theta<theta0:
                Y.append(y)
                return Y
        y.append(P)
        #handles=coating.plotPoint(env, P, handles,array((1,0,0)))
    

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

def nextLine(P0):
    d0 = sqrt(P0[0]**2+P0[1]**2+P0[2]**2)
    R0=1.425
    n = math.ceil((d0-R0)/0.003)
    Rn = R0+n*0.003
    return Rn

def initialPoint():
    rRp = rR[:,0:3]
    P0 = rRp[argmin(rRp[:,1])]
    P0 = nearestSample(P0,rRp)
    Pd=FindNextParallel(P0,1)
    return Pd

def main():
    global handles
    global Rn
    RBF.set_handles(handles,env)
    #env.SetViewer('qtcoin')
    
    try:
        Y = load('trajectory/Y.npz')
        Y = Y['array']
        tempY = []
        for y in Y:
            tempY.append(list(y))   
        Y = tempY
        Pd = Y[-1][-1]
        Rn = math.sqrt(Pd[0]**2+Pd[1]**2+Pd[2]**2)-loopdistance*0.003
        Pd = curvepoint([Pd[0],Pd[1],Pd[2]])
        norm = RBF.df(Pd)
        norm /= sqrt(dot(norm,norm))
        Pd = [Pd[0],Pd[1],Pd[2],norm[0],norm[1],norm[2]]
        #coating.plotPointsArray(env, Y, handles,array((0,1,0)))
    except:
        Y=[[]]
        Pd = initialPoint()
        
    while Rn>1.425:
        Y = drawParallel(Y,Pd)
        savez_compressed('trajectory/'+'Y.npz', array=Y)   
        P0=Y[-1][-1]
        Rn-=loopdistance*0.003

        Pd = curvepoint([P0[0],P0[1],P0[2]])    
        norm = RBF.df(Pd)
        norm /= sqrt(dot(norm,norm))
        Pd = [Pd[0],Pd[1],Pd[2],norm[0],norm[1],norm[2]]
    return Y

if __name__ == '__main__':
    Y = main()
    savez_compressed('trajectory/'+'Y.npz', array=Y)  
