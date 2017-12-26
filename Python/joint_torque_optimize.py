# -*- coding: utf-8 -*-
import coating
from numpy import *
from openravepy import *
from openravepy.misc import SpaceSamplerExtra
import scipy
import random
from scipy.spatial import KDTree

env=Environment()
#env.SetViewer('qtcoin')
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
y = 0.5

# LOAD SAMPLES
indexBlackList = load('indexBlackList_0d_pa24_c12_t30_y05_b1.npz')
indexBlackList  = indexBlackList['array']
indexBlackList = indexBlackList[0]
indexBlueList = load('indexBlueList_0d_pa24_c12_t30_y05_b1.npz')
indexBlueList  = indexBlueList['array']
indexBlueList = indexBlueList[0]
reachableRays = load('reachableRays_0d_pa24_c12_t30_y05_b1.npz')
reachableRays  = reachableRays['array']
AllreachableRays = load('AllreachableRays_0d_pa24_c12_t30_y05_b1.npz')
AllreachableRays  = AllreachableRays['array']
iksolList = load('iksolList_0d_pa24_c12_t30_y05_b1.npz')
iksolList  = iksolList['array']
AlliksolList = load('AlliksolList_0d_pa24_c12_t30_y05_b1.npz')
AlliksolList  = AlliksolList['array']

# 24 DEGREES BASE POSITION
gbase_position=array([ [ -2.14387087e+00,  -3.22000000e+00,  -9.97857391e-01],
                       [ -1.98593217e+00,  -3.22000000e+00,  -6.64834872e-01],
                       [ -1.82799346e+00,  -3.22000000e+00,  -3.31812352e-01],
                       [ -1.67005475e+00,  -3.22000000e+00,   1.21016736e-03],
                       [ -1.51211605e+00,  -3.22000000e+00,   3.34232687e-01],
                       [ -1.35417734e+00,  -3.22000000e+00,   6.67255206e-01],
                       [ -1.19623863e+00,  -3.22000000e+00,   1.00027773e+00]])

# ROBOT POSITION BASE
base = gbase_position[0]
pN = base
pN[1]+=y
normal = [-1,0,0]
pN = numpy.concatenate((pN,normal))
Tn = coating.poseDummy(pN,0)
robot.SetTransform(Tn)

#CHANGE ORIENTATION 
Ti = []
for body in env.GetBodies():
    Ti.append(body.GetTransform())

alpha = 1.0*BladePosition*pi/180
T = numpy.array([[1,0,0,0],[0,cos(alpha),-sin(alpha),0],
                 [0,sin(alpha),cos(alpha),0],[0,0,0,1]])

for ibody in range(0,len(env.GetBodies())):
    if ibody!=5:
        env.GetBodies()[ibody].SetTransform(dot(T,Ti[ibody]))


# ANGLE COMPUTATION
def joint_decoupler(diff,joints):
    for dif in diff:
            dif = dif[joints]
    return diff        

def criteria(diffs,joints):
    diffs = array(diffs)
    sums = []
    for diff in diffs:
        diff = joint_decoupler(diff,joints)
        sums.append(sum(diff,0)*1.0/len(diff))
    minisum = 100
    index = 0
    for i in range(0,len(sums)):    
        if max(sums[i])<minisum:
            minisum = max(sums[i])
            index = i
    return index        

def stopcriteria(T,initialsols,joints):
    index = range(0,len(initialsols))
    sums = [0,0]
    while len(index)>0:
        idxs = T.query_ball_point(reachableRays[index[0],0:3],r=0.05)[1:]
        diff = []
        for idx in idxs:
            diff.append((initialsols[k]-initialsols[idx])**2)
        diff = joint_decoupler(array(diff),joints)
        sums = sums + sum(diff,0)
        index.remove(index[0])
    return sums    

rR = concatenate((reachableRays,AllreachableRays))
T = KDTree(rR[:,0:3])

iI = concatenate((iksolList,AlliksolList))
initialsols = []
for i in iI:
    initialsols.append(i[0])

joints = [3,4]
stopcriteria(T,initialsols,joints)
while True:
    Index = range(0,len(iI))
    while len(Index)>0:
        k = random.choice(Index)
        idxs = T.query_ball_point(reachableRays[k,0:3],r=0.05)[1:]
        possibleSols = iI[k]
        diffs = []
        for sol in possibleSols:
            diff = []
            for idx in idxs:
                diff.append((sol-initialsols[idx])**2)
            diffs.append(diff)
        initialsols[k]=possibleSols[criteria(diffs,joints)]
        Index.remove(k)
    sums = stopcriteria(T,initialsols,joints)
    print(sums)
