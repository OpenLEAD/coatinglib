import coating
from numpy import *
from openravepy import *
from openravepy.misc import SpaceSamplerExtra
import scipy

env=Environment()
env.SetViewer('qtcoin')
env.Load("/home/renan/workspace/coatinglib/Turbina/env_mh12_0_16.xml")
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


pN = numpy.array([ -1.7, -3.3, 0 ])
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

approachrays = load('blade_faro_fast.npz')
approachrays = approachrays['array']
coatedarray = numpy.load('coatedpoints0_faro.npz')
coatedarray  = coatedarray['array']

N = approachrays.shape[0]

# Initial position
BladePositions = [-20,-10,10,35]
RobotPositions = [0,0.1,-0.1,-1.5]
handles=[]
index=3
robottobladedistance = RobotPositions[index]
alpha = 1.0*BladePositions[index]*pi/180;
T = numpy.array([[1,0,0,0],[0,cos(alpha),-sin(alpha),0],
                     [0,sin(alpha),cos(alpha),0],[0,0,0,1]])
i=0
for body in env.GetBodies():
    body.SetTransform(dot(T,Ti[i]))
    i+=1

Ttarget = target.GetTransform()
gapproachrays = c_[dot(approachrays[0:N,0:3],transpose(Ttarget[0:3,0:3]))+tile(Ttarget[0:3,3],(N,1)),dot(approachrays[0:N,3:6],transpose(Ttarget[0:3,0:3]))]

handles.append(env.plot3(points=gapproachrays[:,0:3],pointsize=5,colors=array((1,0,0))))
#handles.append(env.plot3(points=coatedarray[:,0:3],pointsize=5,colors=array((0,0,0))))    
Tn = coating.poseDummy(pN,robottobladedistance)
robot.SetTransform(Tn)
