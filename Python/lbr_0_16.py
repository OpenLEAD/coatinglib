import coating
from numpy import *
from openravepy import *
from openravepy.misc import SpaceSamplerExtra
env=Environment()
env.SetViewer('qtcoin')
env.Load("/home/renan/Documents/EMMA/Turbina/env_LBR_0_16.xml")
robot = env.GetRobots()[0]
target = env.GetBodies()[0]
manip = robot.GetActiveManipulator()


# PARAMETERS
facevector = [0,-1,0]
theta = [0,0,0]
coatingdistance = 0.23 # coating distance
robottobladedistance = 0.12 # robot to blade distance
numberofangles = 8 # degree step
tolerance = 30 # degrees


#pN = numpy.array([-1.412,-2.567,-0.617])# Extremo esquerdo superior

#pN = numpy.array([-0.576,-2.984,0.106])# Extremo direito superior

#pN = numpy.array([-0.4,-3.573,0.26])# Extremo direito inferior

pN = numpy.array([ -1.044, -3.218, 0 ])

normal = [-1,0,0]
pN = numpy.concatenate((pN,normal))

#MAIN
with env:
    ikmodel = databases.inversekinematics.InverseKinematicsModel(robot=robot,iktype=IkParameterization.Type.Transform6D)
    if not ikmodel.load():
        ikmodel.autogenerate()

    approachrays = load('bladepoints16Back.npz')
    approachrays = approachrays['array']
    N = approachrays.shape[0]
    Ttarget = target.GetTransform()
   
    #PLOT BLADE POINTS FOR COATING
    gapproachrays = c_[dot(approachrays[0:N,0:3],transpose(Ttarget[0:3,0:3]))+tile(Ttarget[0:3,3],(N,1)),dot(approachrays[0:N,3:6],transpose(Ttarget[0:3,0:3]))]
    approachgraphs = env.plot3(points=gapproachrays[:,0:3],pointsize=5,colors=array((1,0,0)))


    # Compute Solutions
    reachableRays, iksolList = coating.WorkspaceOnPose(pN, robottobladedistance, gapproachrays,robot,ikmodel,facevector,theta)

    #EXTRA COATING
    AllreachableRays, AlliksolList = coating.AllExtraCoating(gapproachrays,reachableRays,coatingdistance,numberofangles,tolerance,ikmodel,facevector)
    
    #PLOT REACHABLE POINT
    if len(reachableRays)>0:
        grays = c_[reachableRays[:,0:3],reachableRays[:,3:6]]
        raygraphs = env.plot3(points=reachableRays[:,0:3],pointsize=5,colors=array((0,0,0)))

    if len(AllreachableRays)>0:
        grays2 = c_[AllreachableRays[:,0:3],AllreachableRays[:,3:6]]
        raygraphs2 = env.plot3(points=AllreachableRays[:,0:3],pointsize=5,colors=array((0,0,1)))

