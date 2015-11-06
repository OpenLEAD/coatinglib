from numpy import *
from openravepy import *
from openravepy.misc import SpaceSamplerExtra
import scipy

env=Environment()
#env.SetViewer('qtcoin')
env.Load("/home/renan/workspace/coatinglib/Turbina/env_mh12_0_16.xml")
target = env.GetBodies()[0]
target.SetTransform(eye(4))

handles = []
approachrays = load('blade_faro_fast2.npz')
approachrays = approachrays['array']

points = approachrays[approachrays[:,1]<-3]
points = points[points[:,1]>-3.1]
points = points[points[:,2]>-1.5]
points = points[points[:,2]<2]

c=-0.3
beginRow = points[points[:,2]==min(points[:,2])][0]
endRow = points[points[:,2]==max(points[:,2])][0]
beginRow = beginRow[:3]+c*beginRow[3:6]
endRow = endRow[:3]+c*endRow[3:6]
endRow[1] = beginRow[1]
Row=[]
for t in range(101):
    Row.append(beginRow+(endRow-beginRow)*0.01*t)
Row = array(Row)

#handles.append(env.plot3(points=points[:,0:3],pointsize=5,colors=array((0,0,0))))
#handles.append(env.plot3(points=Row,pointsize=5,colors=array((1,0,0))))

savez_compressed('base_position.npz', array=Row)

