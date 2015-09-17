import coating
from numpy import *
from openravepy import *
from openravepy.misc import SpaceSamplerExtra
import time

env=Environment()
env.SetViewer('qtcoin')
env.Load("/home/renan/Documents/EMMA/Turbina/env_SIA20D.xml")
robot = env.GetRobots()[0]
target = env.GetBodies()[0]
manip = robot.GetActiveManipulator()

# PARAMETERS
delta = 0.02
normalanglerange = 0
directiondelta = 0.4
anglerange = pi/6
coatingdistance = 0.23
#robottobladedistance = 1.51
robottobladedistance = 1.2
#pN = numpy.array([-0.00912,-2.72560,-0.07085,0.66103,-0.02598,-0.74991])
#pN = numpy.array([0.005,-3.422,-0.035,0.693,-0.02603,-0.72047])
#pN = numpy.array([0.14924, -2.58874, 0.06677, 0.65749, -0.03902, -0.75245])
pN = numpy.array([-0.03656,-2.98442,-0.08606,0.66103,-0.02598,-0.74991])

initialdistance = 1.5
facevector = [0,1,0]
theta = [0,0,0]

#MAIN
with env:
    ikmodel = databases.inversekinematics.InverseKinematicsModel(robot=robot,
                                                                 iktype=IkParameterization.Type.Transform6D)#,freeindices=[0])
    if not ikmodel.load():
        ikmodel.autogenerate()

    approachrays = load('bladepoints2.npz')
    approachrays = approachrays['array']
    N = approachrays.shape[0]
    Ttarget = target.GetTransform()

    I=0
    M=N-I     
    #PLOT BLADE POINTS FOR COATING
    gapproachrays = c_[dot(approachrays[I:N,0:3],transpose(Ttarget[0:3,0:3]))+tile(Ttarget[0:3,3],(M,1)),dot(approachrays[I:N,3:6],transpose(Ttarget[0:3,0:3]))]
    approachgraphs = env.plot3(points=gapproachrays[:,0:3],pointsize=5,colors=array((1,0,0)))
    
    reachableRays, iksolList = coating.WorkspaceOnPose(pN, robottobladedistance, approachrays,robot,ikmodel,facevector,theta)
    
    #PLOT REACHABLE POINT
    grays = c_[dot(reachableRays[:,0:3],transpose(Ttarget[0:3,0:3]))+tile(Ttarget[0:3,3],(len(reachableRays),1)),dot(reachableRays[:,3:6],transpose(Ttarget[0:3,0:3]))]
    raygraphs = env.plot3(points=reachableRays[:,0:3],pointsize=5,colors=array((0,0,0)))
