import coating
from numpy import *
from openravepy import *
from openravepy.misc import SpaceSamplerExtra
import scipy
import time

env=Environment()
env.SetViewer('qtcoin')
env.Load("/home/renan/workspace/coatinglib/Turbina/env_mh12_0_16.xml")

# CAMERA SETTINGS
Tcamera = numpy.array([[0.53056445,0.0478718,0.84629171,-3.62191391],
       [ 0.25713833,-0.96044621,-0.10687824,-2.01647758],
       [ 0.80770121,0.27431983,-0.52188829, 2.76554561],
       [ 0,0,0,1]])
env.GetViewer().SetCamera(Tcamera)

I = env.GetViewer().GetCameraImage(640,480,Tcamera,[640,640,320,240])
scipy.misc.imsave('/home/renan/Documents/EMMA/Python/framebyframe/distance0/openrave.jpg',I)
