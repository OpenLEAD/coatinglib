import coating
from numpy import array, cross
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

_PRIMARY_RAIL = "primary_rail"
_SECONDARY_RAIL = "secondary_rail"


class Turbine:

    def __init__(self,env,floor_origin,nose_axis):
        self._env = env
        self._primary = next(body for body in bodies if body.GetName()==_PRIMARY_RAIL)
        self._secondary = next(body for body in bodies if body.GetName()==_SECONDARY_RAIL)
        self._floor_origin = floor_origin
        self._floor_axis = array([nose_axis, cross(floor_origin,nose_axis)])
        

    def floor2xyz( v ):
        return v
    
    def xyz2floor( v ):
        return v
