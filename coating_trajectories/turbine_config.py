import ConfigParser
from os import path

class ConfigFileError(Exception):
    pass

class TurbineConfig:
    path = None
    source = None
    
    class Environment:
        load = None
        z_floor_level = None
        primary_safe_margin = None
        secondary_safe_margin = None
        robot_level_difference = None
        blade_angle = None
        rotor_angle = None
        x_max = None
        x_min = None
        y_max = None
        y_min = None
        rail_angle_mean = None
        rail_angle_limit = None
        
    class Coating:
        min_distance = None
        ideal_distance = None
        max_distance = None
        angle_tolerance = None
        coating_speed = None
        parallel_gap = None

    class Model:
        nose_radius = None
        runner_radius = None
        trajectory_step = None

    def __init__(self):
        self.environment = TurbineConfig.Environment()
        self.coating = TurbineConfig.Coating()
        self.model = TurbineConfig.Model()
        
    @classmethod
    def load(cls, config_file, dir_path = ""):
        # Standard python file parsing
        
        parser = ConfigParser.RawConfigParser()
        if parser.read(path.join(dir_path,config_file)) == []:
            raise ConfigFileError("No file named "+path.join(dir_path,config_file)+".")
        
        config = cls()

        config.dir_path = dir_path
        config.source = config_file
        
        try: 
            # environment Section
            config.environment.load = parser.get("environment","load")
            config.environment.z_floor_level = parser.getfloat("environment","z_floor_level")
            config.environment.primary_safe_margin = parser.getfloat("environment","primary_safe_margin")
            config.environment.secondary_safe_margin = parser.getfloat("environment","secondary_safe_margin")
            config.environment.robot_level_difference = parser.getfloat("environment","robot_level_difference")
            config.environment.blade_angle = parser.getfloat("environment","blade_angle")
            config.environment.rotor_angle = parser.getfloat("environment","rotor_angle")
            config.environment.x_max = parser.getfloat("environment","x_max")
            config.environment.x_min = parser.getfloat("environment","x_min")
            config.environment.y_max = parser.getfloat("environment","y_max")
            config.environment.y_min = parser.getfloat("environment","y_min")
            config.environment.rail_angle_mean = parser.getfloat("environment","rail_angle_mean")
            config.environment.rail_angle_limit = parser.getfloat("environment","rail_angle_limit")

            # coating Section
            config.coating.min_distance = parser.getfloat("coating","min_distance")
            config.coating.ideal_distance = parser.getfloat("coating","ideal_distance")
            config.coating.max_distance = parser.getfloat("coating","max_distance")
            config.coating.angle_tolerance = parser.getfloat("coating","angle_tolerance")
            config.coating.coating_speed = parser.getfloat("coating","coating_speed")
            config.coating.parallel_gap = parser.getfloat("coating","parallel_gap")

            # model Section
            config.model.nose_radius = parser.getfloat("model","nose_radius")
            config.model.runner_radius = parser.getfloat("model","runner_radius")
            config.model.trajectory_step = parser.getfloat("model","trajectory_step")
            
        except ConfigParser.NoOptionError as error:
            raise ConfigFileError("Missing "+error.option+" on section "+error.section+".")

        return config
