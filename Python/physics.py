import coating
from openravepy import *
import numpy, time


env=Environment()
env.Load("/home/renan/Documents/EMMA/motoman/motoman_mh122p.xml")
env.SetViewer('qtcoin')

robot = env.GetRobots()[0]
DOF = robot.GetDOF()
manip = robot.GetActiveManipulator()
links = robot.GetLinks()
joints = robot.GetJoints()

# Get poinst for coating
coatedarray = numpy.load('coatedpoints.npz')
coatedarray  = coatedarray['array']

coatedsolarray = numpy.load('ikcoatedpoints.npz')
coatedsolarray  = coatedsolarray['array']

# Robot initial position
pN = numpy.array([ -1, -3.3, 0 ])
normal = [-1,0,0]
pN = numpy.concatenate((pN,normal))
distances = [0,0.1,-0.1,-0.3]
distance = distances[0]
Tn = coating.poseDummy(pN, distance)
robot.SetTransform(Tn)
#RobotPositions = [0,0.1,-0.1,-0.3]
minVelocity = 1.0*40/60 #m/s

# ikmodel
ikmodel = databases.inversekinematics.InverseKinematicsModel(robot=robot,iktype=IkParameterization.Type.Transform6D)
if not ikmodel.load():
    ikmodel.autogenerate()

# CAMERA SETTINGS
Tcamera = numpy.array([[ 0.05777387, -0.06852652,  0.99597505, -4.08520365],
       [-0.32092178, -0.94596499, -0.04646983, -1.95519543],
       [ 0.94534194, -0.31694535, -0.07664371, -0.661735  ],
       [ 0        ,  0        ,  0        ,  1        ]])

env.GetViewer().SetCamera(Tcamera)

##with env:
##    # set a physics engine
##    physics = RaveCreatePhysicsEngine(env,'ode')
##    env.SetPhysicsEngine(physics)
##    physics.SetGravity(numpy.array((0,-9.8,0)))
##    links[0].SetStatic(True)
##    env.StopSimulation()
##    env.StartSimulation(timestep=0.001)

#while True:
#    torques = 100*(numpy.random.rand(DOF)-0.5)
#    for i in range(100):
#        robot.SetJointTorques(torques,True)
#        time.sleep(0.01)
