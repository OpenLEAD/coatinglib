import coating
from numpy import *
from openravepy import *
from openravepy.misc import SpaceSamplerExtra
import time

env=Environment()
env.SetViewer('qtcoin')
env.Load("/home/renan/Documents/EMMA/Turbina/env_LBR.xml")
robot = env.GetRobots()[0]
target = env.GetBodies()[0]
manip = robot.GetActiveManipulator()

# PARAMETERS
delta = 0.02
normalanglerange = 0
directiondelta = 0.4
anglerange = pi/6
coatingdistance = 0.23
#robottobladedistance = 1
#robottobladedistance = 0.64
robottobladedistance = 0
#pN = numpy.array([-0.00912,-2.72560,-0.07085,0.66103,-0.02598,-0.74991])
#pN = numpy.array([0.005,-3.422,-0.035,0.693,-0.02603,-0.72047])
#pN = numpy.array([0.14924, -2.58874, 0.06677, 0.65749, -0.03902, -0.75245])
#pN = numpy.array([-1.051,-2.618,-1.618,1,0,0])
#pN = numpy.array([-0.9583,-2.29458,-2.38312,-0.72054,-0.38453,-0.57702])
pN = numpy.array([-0.87,-2.14,-2.17])
pN = numpy.array([0.076, -3.435, 0.279])
facevector = [0,1,0]
theta = [0,0,0]

normal = [1,0,0]
pN = numpy.concatenate((pN,normal))

#MAIN
with env:
    ikmodel = databases.inversekinematics.InverseKinematicsModel(robot=robot,
                                                                 iktype=IkParameterization.Type.Transform6D)#,freeindices=[0])
    if not ikmodel.load():
        ikmodel.autogenerate()

    approachrays = load('bladepointsl.npz')
    approachrays = approachrays['array']
    N = approachrays.shape[0]
    Ttarget = target.GetTransform()

    approachrays2 = load('bladepointsLeft2.npz')
    approachrays2 = approachrays2['array']
    N2 = approachrays2.shape[0]
   
    #PLOT BLADE POINTS FOR COATING
    gapproachrays = c_[dot(approachrays[0:N,0:3],transpose(Ttarget[0:3,0:3]))+tile(Ttarget[0:3,3],(N,1)),dot(approachrays[0:N,3:6],transpose(Ttarget[0:3,0:3]))]
    approachgraphs = env.plot3(points=gapproachrays[:,0:3],pointsize=5,colors=array((1,0,0)))

    gapproachrays2 = c_[dot(approachrays2[0:N2,0:3],transpose(Ttarget[0:3,0:3]))+tile(Ttarget[0:3,3],(N2,1)),dot(approachrays2[0:N2,3:6],transpose(Ttarget[0:3,0:3]))]
    approachgraphs2 = env.plot3(points=gapproachrays2[:,0:3],pointsize=5,colors=array((1,0,0)))
    
    #reachableRays, iksolList = coating.WorkspaceOnPose(pN, robottobladedistance, approachrays,robot,ikmodel,facevector,theta)
    reachableRays2, iksolList2 = coating.WorkspaceOnPose(pN, robottobladedistance, approachrays2,robot,ikmodel,facevector,theta)
    
    #PLOT REACHABLE POINT
    #if len(reachableRays)>0:
        #grays = c_[dot(reachableRays[:,0:3],transpose(Ttarget[0:3,0:3]))+tile(Ttarget[0:3,3],(len(reachableRays),1)),dot(reachableRays[:,3:6],transpose(Ttarget[0:3,0:3]))]
        #raygraphs = env.plot3(points=reachableRays[:,0:3],pointsize=5,colors=array((0,0,0)))
    if len(reachableRays2)>0:
        grays2 = c_[dot(reachableRays2[:,0:3],transpose(Ttarget[0:3,0:3]))+tile(Ttarget[0:3,3],(len(reachableRays2),1)),dot(reachableRays2[:,3:6],transpose(Ttarget[0:3,0:3]))]
        raygraphs2 = env.plot3(points=reachableRays2[:,0:3],pointsize=5,colors=array((0,0,0)))
