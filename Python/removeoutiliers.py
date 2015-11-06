from numpy import *
from openravepy import *
from openravepy.misc import SpaceSamplerExtra
import scipy

env=Environment()
env.SetViewer('qtcoin')
env.Load("/home/renan/workspace/coatinglib/Turbina/env_mh12_0_16.xml")
target = env.GetBodies()[0]
target.SetTransform(eye(4))

handles = []
approachrays = load('blade_faro_fast2.npz')
approachrays = approachrays['array']

points = approachrays[approachrays[:,1]<-2.5]
points=points[points[:,2]>1.7]

handles.append(env.plot3(points=points[:,0:3],pointsize=5,colors=array((1,0,0))))

savez_compressed('blade_foot2.npz', array=points)
