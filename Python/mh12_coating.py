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

# CAMERA SETTINGS
Tcamera = numpy.array([[ 0.05777387, -0.06852652,  0.99597505, -4.08520365],
       [-0.32092178, -0.94596499, -0.04646983, -1.95519543],
       [ 0.94534194, -0.31694535, -0.07664371, -0.661735  ],
       [ 0        ,  0        ,  0        ,  1        ]])

env.GetViewer().SetCamera(Tcamera)

#MAIN
ikmodel = databases.inversekinematics.InverseKinematicsModel(robot=robot,iktype=IkParameterization.Type.Transform6D)
if not ikmodel.load():
    ikmodel.autogenerate()

approachrays = load('bladepoints16Back.npz')
approachrays = approachrays['array']
N = approachrays.shape[0]

indexlist = zeros((len(approachrays),1),dtype=bool)

# Initial position
handles=[]
a=25

alpha = 1.0*a*pi/180;
T = numpy.array([[1,0,0,0],[0,cos(alpha),-sin(alpha),0],
                 [0,sin(alpha),cos(alpha),0],[0,0,0,1]])

for i in range(0,5):
    env.GetBodies()[i].SetTransform(dot(T,env.GetBodies()[i].GetTransform()))

#PLOT BLADE POINTS FOR COATING
Ttarget = target.GetTransform()
gapproachrays = c_[dot(approachrays[0:N,0:3],transpose(Ttarget[0:3,0:3]))+tile(Ttarget[0:3,3],(N,1)),dot(approachrays[0:N,3:6],transpose(Ttarget[0:3,0:3]))]
handles.append(env.plot3(points=gapproachrays[:,0:3],pointsize=5,colors=array((1,0,0))))

# Compute Solutions
reachableRays, iksolList = coating.WorkspaceOnPose(pN, robottobladedistance, gapproachrays,robot,ikmodel,facevector,theta)
handles.append(env.plot3(points=reachableRays[:,0:3],pointsize=5,colors=array((0,0,0))))

#EXTRA COATING
AllreachableRays, AlliksolList = coating.AllExtraCoating(gapproachrays,reachableRays,coatingdistance,numberofangles,tolerance,ikmodel,facevector)
handles.append(env.plot3(points=AllreachableRays[:,0:3],pointsize=5,colors=array((0,0,1))))

env.GetViewer().SendCommand('SetFiguresInCamera 1')
I = env.GetViewer().GetCameraImage(640,480, Tcamera,[640,640,320,240])
scipy.misc.imsave('/home/renan/Documents/EMMA/Python/coating/'+'a25_d03'+'.jpg',I)
    
