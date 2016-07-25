import coating
from numpy import array, eye, concatenate, zeros
from openravepy import transformLookat
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
floor_origin = array([0,-3.22,0])
handles=[]

class Placement:
    """
        This class provides a mecanism to generate uniform random
        positions (in the sense described by the method) for the rail
        and place both robot and rail system on the designed positions.
    """
    def __init__(self, env, robot, target, primary, secondary,
                 floor_origin = array([0,-3.22,0]), primary_safe_margin = 0.2,
                 secondary_safe_margin = 0.2, robot_level_difference = 0.01):
        self._env = env
        self._bodies = env.GetBodies()
        self._target = target
        self._primary = primary
        self._secondary = secondary
        self._robot = robot
        self._floor_origin = floor_origin
        self._primary_safe_margin = primary_safe_margin
        self._secondary_safe_margin = secondary_safe_margin
        self._robot_level_difference = robot_level_difference

    def rail2xyz(self,(p,s,alpha)):
        #P - primary rail,
        #S - secondary rail,
        #alpha - angle from the perpendicular to the primary rail
        return self._floor_origin + array([p+ s*sin(alpha),
                                           0,
                                           s*cos(alpha)])

    def place_robot(self, (p,s,alpha)):
        Placement = eye(4)
        #  offset HARDCODED 0.02 ref Secondary rail collision
        Placement[0:3,3] = self.rail2xyz((p, s, alpha)) + [0, 0.02, 0]
        R = matrixFromAxisAngle([0, alpha, 0])
        robot.SetTransform(dot(Placement, R))

            
    def place_rail(self,(p,s,alpha)):
        #P - primary rail,
        #S - secondary rail,
        #alpha - angle from the perpendicular to the primary rail
        #ISSUE - p and s cannot be ZERO, but ignoring because random
        #Using camera standard axis for rails
        #Get primary and secondary from bodies
        
        # Primary Rail - offset HARDCODED 0.2
        # primary = next(body for body in bodies if body.GetName()=='primary_rail')
        primary_extent = self._primary.GetLinks()[0].GetGeometries()[0].GetBoxExtents()
        primary_offset = array([0,-primary_extent[1],-primary_extent[2]])+array([0, 0, 0.2])
        primary_offset_transform = eye(4)
        primary_offset_transform[0:3,3] = primary_offset
        
        # Secondary Rail - offset HARDCODED 0.2 and 0.01
        # secondary = next(body for body in bodies if body.GetName()=='secondary_rail')
        secondary_extent = self._secondary.GetLinks()[0].GetGeometries()[0].GetBoxExtents()
        #Resizing
        secondary_extent[2] = abs(s)/2.0 + 0.2 
        env.RemoveKinBody(secondary)
        self._secondary.InitFromBoxes(array([concatenate([zeros(3),secondary_extent])]),True)
        env.AddKinBody(secondary)
        #
        secondary_offset = array([0,-secondary_extent[1],secondary_extent[2]])+array([0, 0.01, -0.2])
        secondary_offset_transform = eye(4)
        secondary_offset_transform[0:3,3] = secondary_offset
        
        # Rails Traonsform and Placement
        primary_transform = transformLookat(floor_origin,floor_origin + array([p,0,0]),[0,1,0])
        primary.SetTransform(dot(primary_transform,primary_offset_transform))
        
        secondary_transform = transformLookat(rail2xyz((p,s,alpha)),floor_origin + array([p,0,0]),[0,1,0])    
        secondary.SetTransform(dot(secondary_transform,secondary_offset_transform))

    def rand_rail(N = 1):
        #P - primary rail,
        #S - secondary rail,
        #alpha - angle from the perpendicular to the primary rail
        Slimit = 2 # sencondary rail limit range
        xmax = 1.5; xmin = -1 # Reasonable distances from the blade
        zmax = Slimit; zmin = -Slimit

        # Randomness
        x,z,alpha = numpy.random.rand(3,N) 
        # Random limits
        x = (xmax - xmin)*x + xmin# xmin <= x < xmax - same for z
        z = (zmax - zmin)*z + zmin
        alphalimit = min(arccos(abs(z/Slimit)),pi/6.0) # min( Physical limits, Reasonable Angles)
        alpha = alphalimit*(2*alpha-1) # - alphalimit <= alpha < alphalimit

        # S/P conversion
        S = -z/cos(alpha)
        P = x - S*sin(alpha)
        
        return (P,S,alpha)
    

# Robot Position

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
env.GetViewer().SetCamera(transformLookat(floor_origin +[0,0.5,0],[6,0.5,0],floor_origin+[6,1,0]))

#MAIN
ikmodel = databases.inversekinematics.InverseKinematicsModel(robot=robot,iktype=IkParameterization.Type.Transform6D)
if not ikmodel.load():
    ikmodel.autogenerate()


# Initial position
alpha = 0.10038996418102937 #offset center blade
coating.RotateBodies(env,alpha)
BladePositions = [0]

RobotPositions = [0.9]
p0=array([1, -3.22, 0, 1, 0, 0])

# Initializing points

approachrays = load('blade_sampling_full/blade_crop_ufast.npz')
approachrays = approachrays['array']
N = approachrays.shape[0]

Ttarget = target.GetTransform()
gapproachrays = c_[dot(approachrays[0:N,0:3],transpose(Ttarget[0:3,0:3]))+tile(Ttarget[0:3,3],(N,1)),dot(approachrays[0:N,3:6],transpose(Ttarget[0:3,0:3]))]

indexlistblack = zeros((len(approachrays),),dtype=bool)
indexlistblue = zeros((len(approachrays),),dtype=bool)


# Randomize and check collisions
bodies = env.GetBodies()
primary = next(body for body in bodies if body.GetName()=='primary_rail')
secondary = next(body for body in bodies if body.GetName()=='secondary_rail')

while True:
    rd = rand_rail()
    place_rail(rd)
    if env.CheckCollision(primary,target):
        print 'Primary rail collision... Redoing'
        continue
    if env.CheckCollision(secondary,target):
        print 'Secondary rail collision... Redoing'
        continue
    
    place_robot(rd)
    collisions = [env.CheckCollision(robot,body) for body in bodies];
    if True in collisions:
        print 'Robot collision with ', bodies[collisions.index(True)].GetName(),'... Redoing'
        continue
    
    break

for pos in BladePositions:

    for robottobladedistance in RobotPositions:
        # Compute Solutions
        #reachableRays, iksolList, indexlist1 = coating.WorkspaceOnPose(p0, robottobladedistance, gapproachrays,robot,ikmodel,facevector,theta)
        reachableRays, iksolList, indexlist1 = coating.Workspace(gapproachrays,robot,ikmodel,facevector,theta)
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


def try_new_rail():
    while True:
        rd = rand_rail()
        place_rail(rd)
        if env.CheckCollision(primary,target):
            print 'Primary rail collision... Redoing'
            continue
        if env.CheckCollision(secondary,target):
            print 'Secondary rail collision... Redoing'
            continue
        
        place_robot(rd)
        collisions = [env.CheckCollision(robot,body) for body in bodies];
        if True in collisions:
            print 'Robot collision with ', bodies[collisions.index(True)].GetName(),'... Redoing'
            continue
    
        break

    Ttarget = target.GetTransform()
    gapproachrays = c_[dot(approachrays[0:N,0:3],transpose(Ttarget[0:3,0:3]))+tile(Ttarget[0:3,3],(N,1)),dot(approachrays[0:N,3:6],transpose(Ttarget[0:3,0:3]))]

    indexlistblack = zeros((len(approachrays),),dtype=bool)
    indexlistblue = zeros((len(approachrays),),dtype=bool)
    # Compute Solutions
    reachableRays, iksolList, indexlist1 = coating.Workspace(gapproachrays,robot,ikmodel,facevector,theta)
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

