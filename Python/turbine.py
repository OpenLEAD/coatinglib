import ConfigParser
from numpy import array, cross, dot
from openravepy import Environment
##from openravepy.misc import SpaceSamplerExtra
##import scipy
##from random import *
##import sys



_PARSE_SECTION = "Environment"
_PARSE_LOAD = "load"
_PARSE_FLOOR_ORIGIN = "floor_origin"
_PARSE_NOSE_AXIS = "nose_axis"

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
                 _RUNNER_AREA: "runner area",
                 _IRIS: "iris"
                 }

    def __init__(self,config_file, viewer = True): #env,floor_origin,nose_axis

        config = ConfigParser.RawConfigParser()
        config.read(config_file)
        load = config.get(_PARSE_SECTION,_PARSE_LOAD)

        floor_origin = config.get(_PARSE_SECTION,_PARSE_FLOOR_ORIGIN)
        self._floor_origin = array([float(c) for c in floor_origin.strip('()[]').split(',')])

        nose_axis = config.get(_PARSE_SECTION,_PARSE_NOSE_AXIS)
        self._nose_axis = array([float(c) for c in nose_axis.strip('()[]').split(',')])
        
        
        self.env = Environment()
        self.env.Load("../Turbina/env_mh12_0_16.xml")
        if viewer:
            self.env.SetViewer('qtcoin')
        self.robot = self.env.GetRobots()[0]
        self.manipulator = self.robot.GetActiveManipulator()
        
        
        bodies = self.env.GetBodies()
        
        self._floor_axis = array([self._nose_axis, cross(self._floor_origin,self._nose_axis)])

        
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
