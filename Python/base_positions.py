#!/usr/bin/env python
# -*- coding: utf-8 -*-
import coating
from numpy import *
from openravepy import *
from openravepy.misc import SpaceSamplerExtra
import scipy
import copy
import time

env=Environment()
env.SetViewer('qtcoin')
env.Load("/home/renan/workspace/coatinglib/Turbina/env_mh12_0_16.xml")
robot = env.GetRobots()[0]
target = env.GetBodies()[0]
manip = robot.GetActiveManipulator()
Ttarget = target.GetTransform()
zmin = -1.24
zmax = 1.24

p1 = array([-1.3, -3.23, 0])
p2 = array([1, -3.23, 0])
def transportRail(p1,p2,row_discretization=10):
    Row=[]
    for t in range(row_discretization+1):
        Row.append(p1+(p2-p1)*(1.0/row_discretization)*t)
    Row = array(Row)
    return Row

approachrays = load('blade_sampling/blade_faro_fast2.npz')
approachrays = approachrays['array']
N = approachrays.shape[0]

target.SetTransform(eye(4))

points = approachrays[approachrays[:,1]<-3]
points = points[points[:,1]>-3.1]
points = points[points[:,2]>-1.5]
points = points[points[:,2]<2]

beginRow = points[points[:,2]==min(points[:,2])][0]
endRow = points[points[:,2]==max(points[:,2])][0]

row_discretization = 10
def BaseRow(beginRow=beginRow,endRow=endRow,c=0.5):
    beginPoint = beginRow[:3]+c*beginRow[3:6]
    endPoint = endRow[:3]+c*endRow[3:6]
    endPoint[1] = beginPoint[1]
    Row=[]
    for t in range(row_discretization+1):
        Row.append(beginPoint+(endPoint-beginPoint)*(1.0/row_discretization)*t)
    Row = array(Row)
    return Row

# PARAMETERS
facevector = [1,0,0]
theta = [0,0,0]
coatingdistance = 0.23 # coating distance
coatingdistancetolerance = 0.01
numberofangles = 8 # degree step
tolerance = 30 # degrees
alpha = 1.0*pi/180; #degree blade step
target.SetTransform(Ttarget)

y = [0,0.5]

# Initial T
Ti = []
for body in env.GetBodies():
    Ti.append(body.GetTransform())

#MAIN
ikmodel = databases.inversekinematics.InverseKinematicsModel(robot=robot,iktype=IkParameterization.Type.Transform6D)
if not ikmodel.load():
    ikmodel.autogenerate()

BladePosition = 0
c = 1.2

base_positions = BaseRow(c=c)
##base_positions = transportRail(p1,p2)
base_positions=array(base_positions)

##gbase_position = c_[dot(base_positions[:,0:3],transpose(Ttarget[0:3,0:3]))+tile(Ttarget[0:3,3],(base_positions.shape[0],1))]
##gbase_position = base_positions

alpha = 1.0*BladePosition*pi/180;
T = numpy.array([[1,0,0,0],[0,cos(alpha),-sin(alpha),0],
                 [0,sin(alpha),cos(alpha),0],[0,0,0,1]])

for ibody in range(0,len(env.GetBodies())):
    if ibody!=5:
        env.GetBodies()[ibody].SetTransform(dot(T,Ti[ibody]))
handles=[]
Ttarget = target.GetTransform()
gapproachrays = c_[dot(approachrays[0:N,0:3],transpose(Ttarget[0:3,0:3]))+tile(Ttarget[0:3,3],(N,1)),dot(approachrays[0:N,3:6],transpose(Ttarget[0:3,0:3]))]
handles.append(env.plot3(points=gapproachrays[:,0:3],pointsize=5,colors=array((1,0,0))))
gbase_position = c_[dot(base_positions[:,0:3],transpose(Ttarget[0:3,0:3]))+tile(Ttarget[0:3,3],(base_positions.shape[0],1))]
gbase_position=gbase_position[gbase_position[:,2]>zmin]
gbase_position=gbase_position[gbase_position[:,2]<zmax]
gbase_position[:,1]=-3.22

# posiçao da base para a pa 24 graus
gbase_position=array([[ -2.14387087e+00,  -3.22000000e+00,  -9.97857391e-01],
       [ -1.98593217e+00,  -3.22000000e+00,  -6.64834872e-01],
       [ -1.82799346e+00,  -3.22000000e+00,  -3.31812352e-01],
       [ -1.67005475e+00,  -3.22000000e+00,   1.21016736e-03],
       [ -1.51211605e+00,  -3.22000000e+00,   3.34232687e-01],
       [ -1.35417734e+00,  -3.22000000e+00,   6.67255206e-01],
       [ -1.19623863e+00,  -3.22000000e+00,   1.00027773e+00]])

#gbase_position = base_positions
indexBlack = zeros(len(gapproachrays),dtype=bool)
indexBlue = zeros(len(gapproachrays),dtype=bool)
indexAll = indexBlack|indexBlue
indexList = array([i for i in range(0,len(gapproachrays))])

indexBlackList = []
indexBlueList = []

def completeIndex(indexAll,real_index,indexlist):
    index = real_index[indexlist]
    for i in range(0,len(index)): indexAll[index[i]]=True
    return indexAll
    


discover = []

for dy in y:
    dy_time = time.time()
    for position in range(0,len(gbase_position)):
        position_time = time.time()
        real_index = indexList[~indexAll]
        pN = copy.copy(gbase_position[position])
        pN[1]+=dy
        normal = [-1,0,0]
        pN = numpy.concatenate((pN,normal))
        reachableRays, iksolList, indexlist1 = coating.WorkspaceOnPose(pN, 0, gapproachrays[~indexAll],robot,ikmodel,facevector,theta,coatingdistancetolerance)
        indexBlackList.append(real_index[indexlist1])
        AllreachableRays, AlliksolList, indexlist2 = coating.AllExtraCoating2(gapproachrays[~indexAll],indexlist1,coatingdistance,numberofangles,tolerance,ikmodel,facevector,coatingdistancetolerance)
        indexBlueList.append(real_index[indexlist2])
        try:
            coatedrays = coating.IndexToPoints(gapproachrays[~indexAll],indexlist1|indexlist2)
        except: break    
        indexAll = completeIndex(indexAll,real_index,indexlist1|indexlist2)
        indexBlack = completeIndex(indexBlack,real_index,indexlist1)
        indexBlue = completeIndex(indexBlue,real_index,indexlist2)
        discover.append([dy,position,len(coatedrays)])
        savetxt(str(24)+'.csv',array(discover),fmt='%.4f',delimiter=',')
        print 'position_time = '+str(time.time()-position_time)
    print 'dy_time = '+str(time.time()-dy_time)          

savez_compressed('full_blade_test/indexBlackList_0d_pa24_c12_t30_y05.npz', array=indexBlackList)
savez_compressed('full_blade_test/indexBlueList_0d_pa24_c12_t30_y05.npz', array=indexBlueList)
