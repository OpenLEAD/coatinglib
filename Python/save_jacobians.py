import coating
from openravepy import *
from numpy import *
from math import *

# Get poinst for coating
velocities = load('linearVel0_HD.npz')
velocities = velocities['array']

thetas = load('thetas0_HD.npz')
thetas = thetas['array']

deltasT = load('deltasT0_HD.npz')
deltasT  = deltasT['array']


env=Environment()
env.Load("/home/renan/Documents/EMMA/Turbina/env_mh12_0_16.xml")
robot = env.GetRobots()[0]
target = env.GetBodies()[0]
manip = robot.GetActiveManipulator()



# PARAMETERS
facevector = [1,0,0]
theta = [0,0,0]
coatingdistance = 0.23 # coating distance
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


NewOmegas, Jacobs = coating.calculateOmegasbyJacobian(robot,ikmodel,manip,thetas,velocities,deltasT)

savez_compressed('NewOmegas0_HD.npz', array=NewOmegas)
savez_compressed('Jacobs0_HD.npz', array=Jacobs)
