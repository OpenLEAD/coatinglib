import coating
from openravepy import *
from numpy import *
from math import *

# Get poinst for coating
FeasiblePoints = load('alltriopoints0_HD.npz')
FeasiblePoints  = FeasiblePoints['array']

FeasibleOmegas = load('NewOmegas0_HD.npz')
FeasibleOmegas  = FeasibleOmegas['array']

FeasibleAlphas = load('NewAlphas0_HD.npz')
FeasibleAlphas  = FeasibleAlphas['array']

FeasibleThetas = load('thetas0_HD.npz')
FeasibleThetas = FeasibleThetas['array']

env=Environment()
env.Load("/home/renan/workspace/coatinglib/Turbina/env_mh12_0_16.xml")
robot = env.GetRobots()[0]
target = env.GetBodies()[0]
manip = robot.GetActiveManipulator()
env.GetPhysicsEngine().SetGravity([0,-9.8,0])


# PARAMETERS
facevector = [1,0,0]
theta = [0,0,0]
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
#BladePositions = [-20]
#BladePositions = [-20,-10]
#BladePositions = [-20,-10,10]
BladePositions = [-20,-10,10,30]

#RobotPositions = [0]
#RobotPositions = [0,0.1]
#RobotPositions = [0,0.1,-0.1]
RobotPositions = [0,0.1,-0.1,-0.3]

robottobladedistance = RobotPositions[0]
alpha = 1.0*BladePositions[0]*pi/180;

T = array([[1,0,0,0],[0,cos(alpha),-sin(alpha),0],
                 [0,sin(alpha),cos(alpha),0],[0,0,0,1]])

i=0
for body in env.GetBodies():
    body.SetTransform(dot(T,Ti[i]))
    i+=1

n = len(FeasibleThetas)
Torques = []
for i in range(0,len(FeasibleThetas)):
    Torque = []
    for j in range(0,len(FeasibleThetas[i])):
        robot.SetDOFValues(FeasibleThetas[i][j][1])
        omega = (FeasibleOmegas[i][j][0]+FeasibleOmegas[i][j][1])/2
        robot.SetDOFVelocities(omega)
        torques = robot.ComputeInverseDynamics(FeasibleAlphas[i][j])
        Torque.append(torques)
    Torques.append(Torque)
    print str(i)+'/'+str(n)
    
savez_compressed('NewTorques0_HD.npz', array=Torques)
