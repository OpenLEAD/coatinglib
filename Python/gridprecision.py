import coating
from numpy import *
from openravepy import *
from openravepy.misc import SpaceSamplerExtra
import scipy
from random import *
import sys

env=Environment()
env.SetViewer('qtcoin')
env.Load("../Turbina/env_mh12_0_16.xml")
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


#pN = numpy.array([-1.412,-2.567,-0.617])# Extremo esquerdo superior

#pN = numpy.array([-0.576,-2.984,0.106])# Extremo direito superior

#pN = numpy.array([-0.4,-3.573,0.26])# Extremo direito inferior

pN = numpy.array([ -1, -3.3, 0 ])

normal = [-1,0,0]
pN = numpy.concatenate((pN,normal))

# CAMERA SETTINGS
Tcamera = numpy.array([[ 0.05777387, -0.06852652,  0.99597505, -4.08520365],
       [-0.32092178, -0.94596499, -0.04646983, -1.95519543],
       [ 0.94534194, -0.31694535, -0.07664371, -0.661735  ],
       [ 0        ,  0        ,  0        ,  1        ]])

env.GetViewer().SetCamera(Tcamera)

# Initial T
Ti = []
for body in env.GetBodies():
    Ti.append(body.GetTransform())

#MAIN
ikmodel = databases.inversekinematics.InverseKinematicsModel(robot=robot,iktype=IkParameterization.Type.Transform6D)
if not ikmodel.load():
    ikmodel.autogenerate()

#approachrays = load('bladepoints16Back.npz')
approachrays = load('blade_sampling/blade_crop_fast.npz')
approachrays = approachrays['array']
N = approachrays.shape[0]

indexlistblack = zeros((len(approachrays),),dtype=bool)
indexlistblue = zeros((len(approachrays),),dtype=bool)

# Initial position
BladePositions = [-0]
#BladePositions = [-20,-10]
#BladePositions = [-20,-10,10]
#BladePositions = [-20,-10,10,30]

RobotPositions = [0.9]
#RobotPositions = [0,0.1]
#RobotPositions = [0,0.1,-0.1]
#RobotPositions = [0,0.1,-0.1,-0.3]
handles=[]

borderxy = array([-0.94670736, -3.22, -1.32629051,-0.68951295,  0, -0.722389])
borderpoint = array([-0.94670736, -2.63243936, -1.32629051])
borderpxy = array([-0.94670736, -3.22, -1.32629051])
bordernormal = array([ -0.68951295,  0.05221131, -0.722389  ])
bordernxy = array([ -0.68951295,  0, -0.722389  ])

index=0
for pos in BladePositions:
    alpha = 1.0*pos*pi/180;
    T = numpy.array([[1,0,0,0],[0,cos(alpha),-sin(alpha),0],
                     [0,sin(alpha),cos(alpha),0],[0,0,0,1]])

    i=0
    for body in env.GetBodies():
        body.SetTransform(dot(T,Ti[i]))
        i+=1
    for robottobladedistance in RobotPositions:
        ##PLOT BLADE POINTS FOR COATING
        Ttarget = target.GetTransform()
        gapproachrays = c_[dot(approachrays[0:N,0:3],transpose(Ttarget[0:3,0:3]))+tile(Ttarget[0:3,3],(N,1)),dot(approachrays[0:N,3:6],transpose(Ttarget[0:3,0:3]))]

        # Compute Solutions
        reachableRays, iksolList, indexlist1 = coating.WorkspaceOnPose(borderxy, robottobladedistance, gapproachrays,robot,ikmodel,facevector,theta)
        #EXTRA COATING
        AllreachableRays, AlliksolList, indexlist2 = coating.AllExtraCoating2(gapproachrays,indexlist1,coatingdistance,numberofangles,tolerance,ikmodel,facevector)

        # Index List
        indexlistblack = indexlistblack|indexlist1
        indexlistblue = indexlistblue|indexlist2
    

approachgraphs = env.plot3(points=gapproachrays[:,0:3],pointsize=5,colors=array((1,0,0)))    
#reachableRays = coating.IndexToPoints(gapproachrays,indexlist)    
#handles.append(env.plot3(points=reachableRays[:,0:3],pointsize=5,colors=array((random(),random(),random()))))    

reachableRays = coating.IndexToPoints(gapproachrays,indexlistblack)
extrareachableRays = coating.IndexToPoints(gapproachrays,indexlistblue)

handles.append(env.plot3(points=extrareachableRays[:,0:3],pointsize=5,colors=array((0,0,1))))    
handles.append(env.plot3(points=reachableRays[:,0:3],pointsize=5,colors=array((0,0,0))))    
 
#reachableRays2 = coating.IndexToPoints(gapproachrays,indexlist2)
#handles.append(env.plot3(points=reachableRays2[:,0:3],pointsize=5,colors=array((0,0,1))))    

