import coating
from openravepy import *
from numpy import *
from math import *

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
BladePositions = [-20]
RobotPositions = [0]
index=0


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

pos = BladePositions[index]
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

# Get poinst for coating
Torques = load('NewTorques0_HD.npz')
Torques = Torques['array']

nearlist = numpy.load('nearPointsByNumberOfPoints0_HD.npz')
nearlist  = nearlist['array']

alltriopoints = load('alltriopoints0_HD.npz')
alltriopoints  = alltriopoints['array']

maxTorques = array([9999,9999,9999,30,30,9999])

x=[];y=[];z=[];nx=[];ny=[];nz=[]
torquefail = []
indexes = []

for i in range(0,len(Torques)):
    io = 1
    for j in range(0,len(Torques[i])):
        if (Torques[i][j]<=maxTorques).all()==True:
            x.append(alltriopoints[i][0][0][0])
            y.append(alltriopoints[i][0][0][1])
            z.append(alltriopoints[i][0][0][2])
            io = 0
            break
    if io:
        nx.append(alltriopoints[i][0][0][0])
        ny.append(alltriopoints[i][0][0][1])
        nz.append(alltriopoints[i][0][0][2])
        torquefail.append(alltriopoints[i][0][0])
        indexes.append(i)

torquefail = array(torquefail)
N = torquefail.shape[0]
Ttarget = target.GetTransform()
#gapproachrays = c_[dot(torquefail[0:N,0:3],transpose(Ttarget[0:3,0:3]))+tile(Ttarget[0:3,3],(N,1)),dot(torquefail[0:N,3:6],transpose(Ttarget[0:3,0:3]))]
approachgraphs = env.plot3(points=torquefail[:,0:3],pointsize=5,colors=array((1,0,0)))

reachableRays, iksolList, indexlist1 = coating.WorkspaceOnPose(pN, robottobladedistance, torquefail,robot,ikmodel,facevector,theta,coatingdistancetolerance)
