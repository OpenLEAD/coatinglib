import ConfigParser
from numpy import array, cross, dot
from openravepy import Environment, databases, IkParameterization
from openravepy.misc import InitOpenRAVELogging 
##from openravepy.misc import SpaceSamplerExtra
##import scipy
##from random import *
##import sys


_PRIMARY_RAIL = "primary_rail"
_SECONDARY_RAIL = "secondary_rail"
_BLADE = "pa"
_RUNNER_AREA = "arocamara"
_IRIS = "iris"


class XMLStructError(Exception):    
    def __init__(self, value):
        Exception.__init__(self,"No body named " + value + " found.")
class ConfigFileError(Exception):
    pass


class Turbine:
    
    std_name = { _PRIMARY_RAIL: "primary rail",
                 _SECONDARY_RAIL: "secondary rail",
                 _BLADE: "blade",
                 _RUNNER_AREA: "runner area",
                 _IRIS: "iris"
                 }

    def __init__(self,config_file, viewer = True):
        self._parse_file(config_file)
        
        
        self.env = Environment()
        InitOpenRAVELogging()
        self.env.Load(self.environment.load)
        if viewer:
            self.env.SetViewer('qtcoin')
        self.robot = self.env.GetRobots()[0]
        self.manipulator = self.robot.GetActiveManipulator()
        
        self.ikmodel = databases.inversekinematics.InverseKinematicsModel(
            robot=self.robot,
            iktype=IkParameterization.Type.Transform6D
            )
        
        if not  self.ikmodel.load():
             self.ikmodel.autogenerate()
            
        self.bodies = self.env.GetBodies()
        
        try:
            self.primary = next(body for body in self.bodies if body.GetName()==_PRIMARY_RAIL)
        except StopIteration:
            raise XMLStructError(_PRIMARY_RAIL)

        try:
            self.secondary = next(body for body in self.bodies if body.GetName()==_SECONDARY_RAIL)
        except StopIteration:
            raise XMLStructError(_SECONDARY_RAIL)

        blade_number = 1
        try:
            self.blades = [next(body for body in self.bodies if body.GetName()==_BLADE+str(blade_number))]
        except StopIteration:
            raise XMLStructError(_BLADE+str(blade_number))

        blade_found = True
        while blade_found:
            blade_number += 1
            try:
                self.blades.append(next(body for body in self.bodies if body.GetName()==_BLADE+str(blade_number)))
            except StopIteration:
                blade_found = False

    def _parse_file(self,config_file):
        
        config = ConfigParser.RawConfigParser()
        if config.read(config_file) == []:
            raise ConfigFileError("No file named "+config_file+".")

        class struct:
            pass
        
        self.environment = struct()
        self.coating = struct()
        self.model = struct()
        
        try: 
            # environment Section
            self.environment.load = config.get("environment","load")
            self.environment.z_floor_level = config.getfloat("environment","z_floor_level")
            self.environment.primary_safe_margin = config.getfloat("environment","primary_safe_margin")
            self.environment.secondary_safe_margin = config.getfloat("environment","secondary_safe_margin")
            self.environment.robot_level_difference = config.getfloat("environment","robot_level_difference")
            self.environment.blade_angle = config.getfloat("environment","blade_angle")
            self.environment.rotor_angle = config.getfloat("environment","rotor_angle")

            # coating Section
            self.coating.min_distance = config.getfloat("coating","min_distance")
            self.coating.ideal_distance = config.getfloat("coating","ideal_distance")
            self.coating.max_distance = config.getfloat("coating","max_distance")
            self.coating.angle_tolerance = config.getfloat("coating","angle_tolerance")
            self.coating.coating_speed = config.getfloat("coating","coating_speed")
            self.coating.parallel_gap = config.getfloat("coating","parallel_gap")

            # model Section
            self.model.nose_radius = config.getfloat("model","nose_radius")
            self.model.runner_radius = config.getfloat("model","runner_radius")
            
        except ConfigParser.NoOptionError as error:
            raise ConfigFileError("Missing "+error.option+" on section "+error.section+".")
