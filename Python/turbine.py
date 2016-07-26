import coating
from numpy import array, cross, dot
from openravepy import Environment
##from openravepy.misc import SpaceSamplerExtra
##import scipy
##from random import *
##import sys

env=Environment()
#env.SetViewer('qtcoin')
env.Load("../Turbina/env_mh12_0_16.xml")
robot = env.GetRobots()[0]
target = env.GetBodies()[0]
manip = robot.GetActiveManipulator()
floor_origin = array([0,-3.22,0])
handles=[]

_PRIMARY_RAIL = "primary_rail"
_SECONDARY_RAIL = "secondary_rail"
_BLADE = "pa"
_RUNNER_AREA = "arocamara"
_IRIS = "iris"


class XMLStructError(Exception):    
    def __init__(self, value):
        Exception.__init__(self,"No body named " + value + " found.")


class Turbine:
    
    std_name = { _PRIMARY_RAIL: "primary rail",
                 _SECONDARY_RAIL: "secondary rail",
                 _BLADE: "blade",
                 _RUNNER_AREA: "runner area"
                 _IRIS: "iris"
                 }

    def __init__(self,env,floor_origin,nose_axis):
        self.env = env
        bodies = env.GetBodies()
            
        self._floor_origin = floor_origin
        self._floor_axis = array([nose_axis, cross(floor_origin,nose_axis)])

        
        try:
            self.primary = next(body for body in bodies if body.GetName()==_PRIMARY_RAIL)
        except StopIteration:
            raise XMLStructError(_PRIMARY_RAIL)

        try:
            self.secondary = next(body for body in bodies if body.GetName()==_SECONDARY_RAIL)
        except StopIteration:
            raise XMLStructError(_SECONDARY_RAIL)

        blade_number = 1
        try:
            self.blades = [next(body for body in bodies if body.GetName()==_BLADE+str(blade_number))]
        except StopIteration:
            raise XMLStructError(_BLADE+str(blade_number))

        blade_found = True
        while blade_found:
            blade_number += 1
            try:
                self.blades.append(next(body for body in bodies if body.GetName()==_BLADE+str(blade_number)))
            except StopIteration:
                blade_found = False


    def getFloorOrigin(self):
        return array(self._floor_origin)

    def getFloorAxis(self):
        return array(self._floor_axis)
