from numpy import *
import coating
from openravepy import *
from openravepy.misc import SpaceSamplerExtra
import robot_path

env=Environment()
env.Load("../Turbina/env_mh12_0_16.xml")
robot = env.GetRobots()[0]
target = env.GetBodies()[0]
manip = robot.GetActiveManipulator()
handles=[]
env.SetViewer('qtcoin')
pN=[-1.51211605e+00,-3.22000000e+00,3.34232687e-01]
normal = [-1,0,0]
pN = concatenate((pN,normal))

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
        y = y[y[:,0]<-0.8]
        croppedY.append(y)

    coating.plotPointsArray(env, croppedY, handles,array((0,0,1)))
    savez_compressed('trajectory/'+'croppedY.npz', array=croppedY)
    
Q = robot_path.main(pN, 0, croppedY)
