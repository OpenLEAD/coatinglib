import coating
from numpy import *
from openravepy import *
from openravepy.misc import SpaceSamplerExtra
import scipy
from random import *
import sys

env=Environment()
env.SetViewer('qtcoin')
env.Load("../Turbina/env_mh12_0_16.xml")
robot = env.GetRobots()[0]
target = env.GetBodies()[0]
manip = robot.GetActiveManipulator()
floor = array([0,-3.22,0])
handles=[]

def rail2cart((p,s,alpha)):
    return array([p+s*sin(alpha),0,-s*cos(alpha)])+floor

def place_rail((p,s,alpha)):
    #Get primary and secondary from bodies
    bodies = env.GetBodies()
    
    # Primary Rail - offset HARDCODED
    primary = next(body for body in bodies if body.GetName()=='primary_rail')
    primary_extent = primary.GetLinks()[0].GetGeometries()[0].GetBoxExtents()
    primary_offset = array([-primary_extent[0],0,-primary_extent[2]])+array([0.2 0 0])
    primary_offset_transform = eye(4)
    primary_offset_transform[0:3,3] = primary_offset
    
    # Secondary Rail - offset HARDCODED
    secondary = next(body for body in bodies if body.GetName()=='secondary_rail')
    secondary_extent = secondary.GetLinks()[0].GetGeometries()[0].GetBoxExtents()
    #Resizing
    secondary_extent[0] = (abs(s)+0.2)/2 
    env.RemoveKinBody(secondary)
    secondary.InitFromBoxes(numpy.array([concatenate([zeros(3),secondary_extent])]),True)
    env.AddKinBody(secondary)
    #
    secondary_offset = array([-secondary_extent[0],0,-secondary_extent[2]])+array([0.2 0 0])
    secondary_offset_transform = eye(4)
    secondary_offset_transform[0:3,3] = secondary_offset

    #TODO transfromlookto para colocar os dois rails

def rand_rail(N = 1): #P - primary rail, S - secondary rail, alpha - angle from the perpendicular to the primary rail
    Slimit = 1 # sencondary rail limit range
    xmax = 1.5; xmin = -1
    zmax = Slimit; zmin = -Slimit

    # Randomness
    x,z,alpha = numpy.random.rand(3,N) 
    # Random limits
    x = (xmax - xmin)*x + xmin# xmin <= x < xmax - same for z
    z = (zmax - zmin)*z + zmin
    alphalimit = arccos(abs(z/Slimit))
    alpha = alphalimit*(2*alpha-1) # - alphalimit <= alpha < alphalimit

    # S/P conversion
    S = -z/cos(alpha)
    P = x - S*sin(alpha)
    
    return (P,S,alpha)

# Robot Position
#prime_rail = 

# PARAMETERS
facevector = [1,0,0]
theta = [0,0,0]
coatingdistance = 0.23 # coating distance
robottobladedistance = -0.3 # robot to blade distance
numberofangles = 8 # degree step
tolerance = 20 # degrees

pN = numpy.array([ -1, -3.3, 0 ])

normal = [-1,0,0]
pN = numpy.concatenate((pN,normal))

# CAMERA SETTINGS
env.GetViewer().SetCamera(transformLookat(floor +[0,0.5,0],[6,0.5,0],[0,-1,0]))

#MAIN
ikmodel = databases.inversekinematics.InverseKinematicsModel(robot=robot,iktype=IkParameterization.Type.Transform6D)
if not ikmodel.load():
    ikmodel.autogenerate()

approachrays = load('blade_sampling_full/blade_crop_ufast.npz')
approachrays = approachrays['array']
N = approachrays.shape[0]

indexlistblack = zeros((len(approachrays),),dtype=bool)
indexlistblue = zeros((len(approachrays),),dtype=bool)

# Initial position
alpha = 0.10038996418102937 #offset center blade
coating.RotateBodies(env,alpha)
BladePositions = [0]

RobotPositions = [0.9]
p0=array([1, -3.22, 0, 1, 0, 0])

for pos in BladePositions:

    for robottobladedistance in RobotPositions:
        ##PLOT BLADE POINTS FOR COATING
        Ttarget = target.GetTransform()
        gapproachrays = c_[dot(approachrays[0:N,0:3],transpose(Ttarget[0:3,0:3]))+tile(Ttarget[0:3,3],(N,1)),dot(approachrays[0:N,3:6],transpose(Ttarget[0:3,0:3]))]

        # Compute Solutions
        reachableRays, iksolList, indexlist1 = coating.WorkspaceOnPose(p0, robottobladedistance, gapproachrays,robot,ikmodel,facevector,theta)
        #EXTRA COATING
        AllreachableRays, AlliksolList, indexlist2 = coating.AllExtraCoating2(gapproachrays,indexlist1,coatingdistance,numberofangles,tolerance,ikmodel,facevector)

        # Index List
        indexlistblack = indexlistblack|indexlist1
        indexlistblue = indexlistblue|indexlist2
    

approachgraphs = env.plot3(points=gapproachrays[:,0:3],pointsize=5,colors=array((1,0,0)))    
#reachableRays = coating.IndexToPoints(gapproachrays,indexlist)    
#handles.append(env.plot3(points=reachableRays[:,0:3],pointsize=5,colors=array((random(),random(),random()))))    

reachableRays = coating.IndexToPoints(gapproachrays,indexlistblack)
extrareachableRays = coating.IndexToPoints(gapproachrays,indexlistblue)

handles.append(env.plot3(points=extrareachableRays[:,0:3],pointsize=5,colors=array((0,0,1))))    
handles.append(env.plot3(points=reachableRays[:,0:3],pointsize=5,colors=array((0,0,0))))    
     
#reachableRays2 = coating.IndexToPoints(gapproachrays,indexlist2)
#handles.append(env.plot3(points=reachableRays2[:,0:3],pointsize=5,colors=array((0,0,1))))    

