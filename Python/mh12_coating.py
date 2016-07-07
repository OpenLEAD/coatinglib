# -*- coding: utf-8 -*-
import coating
from numpy import *
from openravepy import *
from openravepy.misc import SpaceSamplerExtra
import scipy
import copy

env=Environment()
env.SetViewer('qtcoin')
env.Load("../Turbina/env_mh12_0_16.xml")
robot = env.GetRobots()[0]
target = env.GetBodies()[0]
manip = robot.GetActiveManipulator()
Ttarget = target.GetTransform()


# PARAMETERS
facevector = [1,0,0]
theta = [0,0,0]
coatingdistance = 0.23 # coating distance
coatingdistancetolerance = 0.01
numberofangles = 8 # degree step
tolerance = 30 # degrees
alpha = 1.0*pi/180; #degree blade step
BladePosition = 0
y = 0

approachrays = load('blade_sampling/blade_crop_fast.npz')
approachrays = approachrays['array']
N = approachrays.shape[0]

#MAIN
ikmodel = databases.inversekinematics.InverseKinematicsModel(robot=robot,iktype=IkParameterization.Type.Transform6D)
if not ikmodel.load():
    ikmodel.autogenerate()

# Initial T
Ti = []
for body in env.GetBodies():
    Ti.append(body.GetTransform())

alpha = 1.0*BladePosition*pi/180;
T = array([[1,0,0,0],[0,cos(alpha),-sin(alpha),0],
                 [0,sin(alpha),cos(alpha),0],[0,0,0,1]])

for ibody in range(0,len(env.GetBodies())):
    if ibody!=5:
        env.GetBodies()[ibody].SetTransform(dot(T,Ti[ibody]))

# Change orientation of samples
Ttarget = target.GetTransform()
gapproachrays = c_[dot(approachrays[0:N,0:3],transpose(Ttarget[0:3,0:3]))+tile(Ttarget[0:3,3],(N,1)),dot(approachrays[0:N,3:6],transpose(Ttarget[0:3,0:3]))]

# posiçao da base para a pa 24 graus
gbase_position=array([[ -2.14387087e+00,  -3.22000000e+00,  -9.97857391e-01],
       [ -1.98593217e+00,  -3.22000000e+00,  -6.64834872e-01],
       [ -1.82799346e+00,  -3.22000000e+00,  -3.31812352e-01],
       [ -1.67005475e+00,  -3.22000000e+00,   1.21016736e-03],
       [ -1.51211605e+00,  -3.22000000e+00,   3.34232687e-01],
       [ -1.35417734e+00,  -3.22000000e+00,   6.67255206e-01],
       [ -1.19623863e+00,  -3.22000000e+00,   1.00027773e+00]])

borderpoint = array([-0.94670736, -2.63243936, -1.32629051])
borderpxy = array([-0.94670736, -3.22, -1.32629051])
bordernormal = array([ -0.68951295,  0.05221131, -0.722389  ])
bordernxy = array([ -0.68951295,  0, -0.722389  ])

ptest = borderpxy+bordernxy

# Saved Variables
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

#PLOT BLADE POINTS FOR COATING
handles=[]
handles.append(env.plot3(points=gapproachrays[:,0:3],pointsize=5,colors=array((1,0,0))))

# Compute Solutions
base = gbase_position[0] # ptest

real_index = indexList[~indexAll]
pN = copy.copy(base)
pN[1]+=y
normal = [-1,0,0]
pN = concatenate((pN,normal))
reachableRays, iksolList, indexlist1 = coating.WorkspaceOnPose(pN, 0, gapproachrays[~indexAll],robot,ikmodel,facevector,theta,coatingdistancetolerance)
indexBlackList.append(real_index[indexlist1])
AllreachableRays, AlliksolList, indexlist2 = coating.AllExtraCoating2(gapproachrays[~indexAll],indexlist1,coatingdistance,numberofangles,tolerance,ikmodel,facevector,coatingdistancetolerance)
indexBlueList.append(real_index[indexlist2])
try:
    coatedrays = coating.IndexToPoints(gapproachrays[~indexAll],indexlist1|indexlist2)
except: None    
indexAll = completeIndex(indexAll,real_index,indexlist1|indexlist2)
indexBlack = completeIndex(indexBlack,real_index,indexlist1)
indexBlue = completeIndex(indexBlue,real_index,indexlist2)          

#savez_compressed('coated_points/pos4_black.npz', array=reachableRays)
#savez_compressed('coated_points/pos4_blue.npz', array=AllreachableRays)

#PLOT

handles.append(env.plot3(points=reachableRays[:,0:3],pointsize=5,colors=array((0,0,0))))
handles.append(env.plot3(points=AllreachableRays[:,0:3],pointsize=5,colors=array((0,0,1))))   
