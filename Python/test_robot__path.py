from numpy import *
import coating
from openravepy import *
from openravepy.misc import SpaceSamplerExtra
import robot_path

pN=[-1.51211605e+00,-3.22000000e+00,3.34232687e-01]
normal = [-1,0,0]
pN = concatenate((pN,normal))
robot_path.env.SetViewer('qtcoin')

try:
    croppedY = load('trajectory/croppedY.npz')
    croppedY = croppedY['array']
except:
    Y = load('trajectory/Y2.npz')
    Y = Y['array']
    croppedY = []
    for y in Y:
        y = array(y)
        y = y[y[:,2]<1.2]
        y = y[y[:,2]>-0.8]
        y = y[y[:,1]>-3]
        y = y[y[:,1]<-1.2]
        y = y[y[:,0]<-0.5]
        croppedY.append(y)

    coating.plotPointsArray(env, croppedY, robot_path.handles,array((0,0,1)))
    savez_compressed('trajectory/'+'croppedY.npz', array=croppedY)

try:
    Q = load('trajectory/Q.npz')
    Q = Q['array']
    croppedTraj = load('trajectory/croppedTraj.npz')
    croppedTraj = croppedTraj['array']
except:  None  
    #Q, croppedTraj = robot_path.main(pN, 0, croppedY)

robot_path.solution(pN, croppedTraj[0][0])
for q in Q:
    coating.robotPath2(q, 0.01,robot_path.robot,robot_path.ikmodel,
                       robot_path.env, robot_path.handles, array((1,0,0)))
