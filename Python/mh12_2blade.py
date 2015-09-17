import coating
from numpy import *
from openravepy import *
from openravepy.misc import SpaceSamplerExtra
import scipy
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
robottobladedistance = 0 # robot to blade distance
numberofangles = 8 # degree step
tolerance = 20 # degrees


#pN = numpy.array([-1.412,-2.567,-0.617])# Extremo esquerdo superior

#pN = numpy.array([-0.576,-2.984,0.106])# Extremo direito superior

#pN = numpy.array([-0.4,-3.573,0.26])# Extremo direito inferior

pN = numpy.array([ -1.044, -3.218, 0 ])

normal = [-1,0,0]
pN = numpy.concatenate((pN,normal))

# CAMERA SETTINGS
Tcamera = numpy.array([[0.53056445,0.0478718,0.84629171,-3.62191391],
       [ 0.25713833,-0.96044621,-0.10687824,-2.01647758],
       [ 0.80770121,0.27431983,-0.52188829, 2.76554561],
       [ 0,0,0,1]])
env.GetViewer().SetCamera(Tcamera)

#MAIN
ikmodel = databases.inversekinematics.InverseKinematicsModel(robot=robot,iktype=IkParameterization.Type.Transform6D)
if not ikmodel.load():
    ikmodel.autogenerate()

approachrays = load('bladepoints16Back.npz')
approachrays = approachrays['array']
N = approachrays.shape[0]

# Initial position
a=15
alpha = 1.0*a*pi/180;
T = numpy.array([[1,0,0,0],[0,cos(alpha),-sin(alpha),0],
                 [0,sin(alpha),cos(alpha),0],[0,0,0,1]])

for i in range(0,5):
    env.GetBodies()[i].SetTransform(dot(T,env.GetBodies()[i].GetTransform()))

    
# COATING LOOP
a=40
while(a):
    alpha = 1.0*pi/180;
    T = numpy.array([[1,0,0,0],[0,cos(alpha),-sin(alpha),0],
                     [0,sin(alpha),cos(alpha),0],[0,0,0,1]])

    for i in range(0,5):
        env.GetBodies()[i].SetTransform(dot(T,env.GetBodies()[i].GetTransform()))


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

    env.GetViewer().SendCommand('SetFiguresInCamera 1')
    I = env.GetViewer().GetCameraImage(640,480, Tcamera,[640,640,320,240])
    scipy.misc.imsave('/home/renan/Documents/EMMA/Python/framebyframe/distance0/'+str(a-10)+'.jpg',I)
    print(abs(a-40))
    a-=1
    
    
    

