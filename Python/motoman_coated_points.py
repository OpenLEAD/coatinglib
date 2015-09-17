import coating
from numpy import *
from openravepy import *
from openravepy.misc import SpaceSamplerExtra
import scipy
from random import *

env=Environment()
env.SetViewer('qtcoin')
env.Load("/home/renan/Documents/EMMA/Turbina/env_mh12_0_16.xml")
robot = env.GetRobots()[0]
target = env.GetBodies()[0]
manip = robot.GetActiveManipulator()



# PARAMETERS
facevector = [1,0,0]
theta = [0,0,0]
coatingdistance = 0.23 # coating distance
robottobladedistance = -0.3 # robot to blade distance
numberofangles = 8 # degree step
tolerance = 20 # degrees
alpha = 1.0*pi/180; #degree blade step

pN = numpy.array([ -1, -3.3, 0 ])
normal = [-1,0,0]
pN = numpy.concatenate((pN,normal))

# Initial T
Ti = []
for body in env.GetBodies():
    Ti.append(body.GetTransform())

#MAIN
ikmodel = databases.inversekinematics.InverseKinematicsModel(robot=robot,iktype=IkParameterization.Type.Transform6D)
if not ikmodel.load():
    ikmodel.autogenerate()

approachrays = load('bladepoints16Back.npz')
approachrays = approachrays['array']
N = approachrays.shape[0]

coatedpointsindex = []
coatedpoints = []
iksolution = []
fulliksollistarray = []

# Positions
BladePositions = [-20,-10,10,30]
RobotPositions = [0,0.1,-0.1,-0.3]
handles=[]

index=0
for pos in BladePositions:
    robottobladedistance = RobotPositions[index]
    alpha = 1.0*pos*pi/180;
    T = numpy.array([[1,0,0,0],[0,cos(alpha),-sin(alpha),0],
                     [0,sin(alpha),cos(alpha),0],[0,0,0,1]])

    i=0
    for body in env.GetBodies():
        body.SetTransform(dot(T,Ti[i]))
        i+=1

    #PLOT BLADE POINTS FOR COATING
    Ttarget = target.GetTransform()
    gapproachrays = c_[dot(approachrays[0:N,0:3],transpose(Ttarget[0:3,0:3]))+tile(Ttarget[0:3,3],(N,1)),dot(approachrays[0:N,3:6],transpose(Ttarget[0:3,0:3]))]

    # Compute Solutions
    reachableRays, iksolList, indexlist1 = coating.WorkspaceOnPose(pN, robottobladedistance, gapproachrays,robot,ikmodel,facevector,theta)
    #EXTRA COATING
    AllreachableRays, AlliksolList, indexlist2 = coating.AllExtraCoating2(gapproachrays,indexlist1,coatingdistance,numberofangles,tolerance,ikmodel,facevector)

    # Index List
    fulliksollist = coating.iksolsSort(indexlist1,indexlist2,iksolList,AlliksolList,ikmodel)
    indexlist = indexlist1|indexlist2
    coatedpointsindex.append(indexlist)
    coated = coating.IndexToPoints(gapproachrays,indexlist)
    coatedpoints.append(coated)
    fulliksollistarray.append(fulliksollist)
    index+=1
    
   
savez_compressed('coatedpoints.npz', array=coatedpoints)
savez_compressed('ikcoatedpoints.npz', array=fulliksollistarray)
