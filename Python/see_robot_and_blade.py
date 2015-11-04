import coating
from numpy import *
from openravepy import *
from openravepy.misc import SpaceSamplerExtra
import scipy
import time

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
coatingdistancetolerance = 0.01
robottobladedistance = -0.3 # robot to blade distance
numberofangles = 8 # degree step
tolerance = 20 # degrees
alpha = 1.0*pi/180; #degree blade step


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

BladePositions = [-20,-10,10,30]
RobotPositions = [0,0.1,-0.1,-0.3]
handles=[]

index=2

robottobladedistance = RobotPositions[index]
pos = BladePositions[index]
alpha = 1.0*pos*pi/180;
T = numpy.array([[1,0,0,0],[0,cos(alpha),-sin(alpha),0],
                 [0,sin(alpha),cos(alpha),0],[0,0,0,1]])

i=0
for body in env.GetBodies():
    body.SetTransform(dot(T,Ti[i]))
    i+=1

Tn = coating.poseDummy(pN,robottobladedistance)
robot.SetTransform(Tn)    
