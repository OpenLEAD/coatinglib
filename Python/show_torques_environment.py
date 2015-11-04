import coating
from openravepy import *
from numpy import *
from math import *

# Get poinst for coating
torquefail = load('torquefail0_HD.npz')
torquefail = torquefail['array']

torquetrue = load('torquetrue0_HD.npz')
torquetrue = torquetrue['array']

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
numberofangles = 8 # degree step
tolerance = 20 # degrees
alpha = 1.0*pi/180; #degree blade step

pN = array([ -1, -3.3, 0 ])
normal = [-1,0,0]
pN = concatenate((pN,normal))

# Initial T
Ti = []
for body in env.GetBodies():
    Ti.append(body.GetTransform())

#MAIN
ikmodel = databases.inversekinematics.InverseKinematicsModel(robot=robot,iktype=IkParameterization.Type.Transform6D)
if not ikmodel.load():
    ikmodel.autogenerate()


# Initial position
BladePositions = [-20,-10,10,30]
RobotPositions = [0,0.1,-0.1,-0.3]

robottobladedistance = RobotPositions[0]
alpha = 1.0*BladePositions[0]*pi/180;

T = array([[1,0,0,0],[0,cos(alpha),-sin(alpha),0],
                 [0,sin(alpha),cos(alpha),0],[0,0,0,1]])

i=0
for body in env.GetBodies():
    body.SetTransform(dot(T,Ti[i]))
    i+=1

# Compute Solutions
reachableRays, iksolList, indexlist1 = coating.WorkspaceOnPose(pN, robottobladedistance, torquetrue,robot,ikmodel,facevector,theta,coatingdistancetolerance)
#EXTRA COATING
AllreachableRays, AlliksolList, indexlist2 = coating.AllExtraCoating2(torquetrue,indexlist1,coatingdistance,numberofangles,tolerance,ikmodel,facevector,coatingdistancetolerance)

iksollisttrue = coating.iksolsSort(indexlist1,indexlist2,iksolList,AlliksolList)

# Compute Solutions
reachableRays, iksolList, indexlist1 = coating.WorkspaceOnPose(pN, robottobladedistance, torquefail,robot,ikmodel,facevector,theta,coatingdistancetolerance)
#EXTRA COATING
AllreachableRays, AlliksolList, indexlist2 = coating.AllExtraCoating2(torquefail,indexlist1,coatingdistance,numberofangles,tolerance,ikmodel,facevector,coatingdistancetolerance)

iksollistfail = coating.iksolsSort(indexlist1,indexlist2,iksolList,AlliksolList)

handles=[]
handles.append(env.plot3(points=torquetrue[:,0:3],pointsize=5,colors=array((0,1,0))))
handles.append(env.plot3(points=torquefail[:,0:3],pointsize=10,colors=array((1,0,0))))
