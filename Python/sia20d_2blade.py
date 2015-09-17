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
pN = numpy.array([-0.9,-2.3,-2.3])
#pN2 = numpy.array([-0.1,-2.7,-1.6])
pN2 = numpy.array([0.2,-2.9,-1.9])
#pN3 = numpy.array([0.9,-2.7,-0.9])
pN3 = numpy.array([0.9,-3.3,-0.9])
pN4 = numpy.array([-0.81,-1.33,-1.5])

#normal = [1,0,0]
#normal = [-0.63,-0.43,0.63]
#normal = [-2824.69483731,-3561.77916006,1205.28249254 ]
normal = [-0.72,-0.39,0.56]
normal = normal/sqrt(dot(normal,normal))

normal2 = [-0.7,-0.37,0.6]
normal2 = normal2/sqrt(dot(normal2,normal2))

normal3 = [-0.7,-0.33,0.64]
normal3 = normal3/sqrt(dot(normal3,normal3))

pN = numpy.concatenate((pN,normal))
pN2 = numpy.concatenate((pN2,normal2))
pN3 = numpy.concatenate((pN3,normal3))
pN4 = numpy.concatenate((pN4,normal))

                  
                  
facevector = [0,1,0]
theta = [0,0,0]

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

    robottobladedistance = 1.2
    #reachableRays, iksolList = coating.WorkspaceOnPose(pN, robottobladedistance, approachrays,robot,ikmodel,facevector,theta)
    reachableRays2, iksolList2 = coating.WorkspaceOnPose(pN, robottobladedistance, approachrays2,robot,ikmodel,facevector,theta)

    robottobladedistance = 1.45
    #reachableRays3, iksolList3 = coating.WorkspaceOnPose(pN2, robottobladedistance, approachrays,robot,ikmodel,facevector,theta)
    reachableRays4, iksolList4 = coating.WorkspaceOnPose(pN2, robottobladedistance, approachrays2,robot,ikmodel,facevector,theta)

    robottobladedistance = 1.4
    #reachableRays5, iksolList5 = coating.WorkspaceOnPose(pN3, robottobladedistance, approachrays,robot,ikmodel,facevector,theta)
    reachableRays6, iksolList6 = coating.WorkspaceOnPose(pN3, robottobladedistance, approachrays2,robot,ikmodel,facevector,theta)

    robottobladedistance = 1.4
    #reachableRays7, iksolList7 = coating.WorkspaceOnPose(pN4, robottobladedistance, approachrays,robot,ikmodel,facevector,theta)
    #reachableRays8, iksolList8 = coating.WorkspaceOnPose(pN4, robottobladedistance, approachrays2,robot,ikmodel,facevector,theta)

    
    #PLOT REACHABLE POINT
    #if len(reachableRays)>0:
        #grays = c_[dot(reachableRays[:,0:3],transpose(Ttarget[0:3,0:3]))+tile(Ttarget[0:3,3],(len(reachableRays),1)),dot(reachableRays[:,3:6],transpose(Ttarget[0:3,0:3]))]
        #raygraphs = env.plot3(points=reachableRays[:,0:3],pointsize=5,colors=array((0,0,0)))
    if len(reachableRays2)>0:
        grays2 = c_[dot(reachableRays2[:,0:3],transpose(Ttarget[0:3,0:3]))+tile(Ttarget[0:3,3],(len(reachableRays2),1)),dot(reachableRays2[:,3:6],transpose(Ttarget[0:3,0:3]))]
        raygraphs2 = env.plot3(points=reachableRays2[:,0:3],pointsize=5,colors=array((0,0,0)))

    #if len(reachableRays3)>0:
        #grays3 = c_[dot(reachableRays3[:,0:3],transpose(Ttarget[0:3,0:3]))+tile(Ttarget[0:3,3],(len(reachableRays3),1)),dot(reachableRays3[:,3:6],transpose(Ttarget[0:3,0:3]))]
        #raygraphs3 = env.plot3(points=reachableRays3[:,0:3],pointsize=5,colors=array((0,1,0)))
    if len(reachableRays4)>0:
        grays4 = c_[dot(reachableRays4[:,0:3],transpose(Ttarget[0:3,0:3]))+tile(Ttarget[0:3,3],(len(reachableRays4),1)),dot(reachableRays4[:,3:6],transpose(Ttarget[0:3,0:3]))]
        raygraphs4 = env.plot3(points=reachableRays4[:,0:3],pointsize=5,colors=array((0,1,0)))

    #if len(reachableRays5)>0:
        #grays5 = c_[dot(reachableRays5[:,0:3],transpose(Ttarget[0:3,0:3]))+tile(Ttarget[0:3,3],(len(reachableRays5),1)),dot(reachableRays5[:,3:6],transpose(Ttarget[0:3,0:3]))]
        #raygraphs5 = env.plot3(points=reachableRays5[:,0:3],pointsize=5,colors=array((0,0,1)))
    if len(reachableRays6)>0:
        grays6 = c_[dot(reachableRays6[:,0:3],transpose(Ttarget[0:3,0:3]))+tile(Ttarget[0:3,3],(len(reachableRays6),1)),dot(reachableRays6[:,3:6],transpose(Ttarget[0:3,0:3]))]
        raygraphs6 = env.plot3(points=reachableRays6[:,0:3],pointsize=5,colors=array((0,0,1)))

    #if len(reachableRays7)>0:
        #grays7 = c_[dot(reachableRays7[:,0:3],transpose(Ttarget[0:3,0:3]))+tile(Ttarget[0:3,3],(len(reachableRays7),1)),dot(reachableRays7[:,3:6],transpose(Ttarget[0:3,0:3]))]
        #raygraphs7 = env.plot3(points=reachableRays7[:,0:3],pointsize=5,colors=array((1,0,1)))
    #if len(reachableRays8)>0:
        #grays8 = c_[dot(reachableRays8[:,0:3],transpose(Ttarget[0:3,0:3]))+tile(Ttarget[0:3,3],(len(reachableRays8),1)),dot(reachableRays8[:,3:6],transpose(Ttarget[0:3,0:3]))]
        #raygraphs8 = env.plot3(points=reachableRays8[:,0:3],pointsize=5,colors=array((1,0,1)))  
