from numpy import array, dot, eye, concatenate, zeros
from openravepy import transformLookat, matrixFromAxisAngle, Environment, databases, IkParameterization
from openravepy.misc import InitOpenRAVELogging
from os import path


_PRIMARY_RAIL = "primary_rail"
_SECONDARY_RAIL = "secondary_rail"
_BLADE = "pa"
_RUNNER_AREA = "arocamara"
_IRIS = "iris"
_ROTOR = "rotor"


class XMLStructError(Exception):    
    def __init__(self, value):
        Exception.__init__(self,"No body named " + value + " found.")



class Turbine:
    """
    This class provides a configuration file loading for the whole turbine environment,
    plus manipulation of the active elements of the turbine (robot/rail)
    and their interactiong (collision)
    """

    env = None
    config = None

    # std names for printing propouses
    std_name = { _PRIMARY_RAIL: "primary rail",
                 _SECONDARY_RAIL: "secondary rail",
                 _BLADE: "blade",
                 _RUNNER_AREA: "runner area",
                 _IRIS: "iris",
                 _ROTOR: "rotor"
                 }

    def __init__(self,config):
        self.config = config
        
        # Create OpenRave environment
        self.env = Environment()
        InitOpenRAVELogging()
        self.env.Load(path.join(config.dir_path,config.environment.load))
        self.robot = self.env.GetRobots()[0]
        self.manipulator = self.robot.GetActiveManipulator()
        
        self.ikmodel = databases.inversekinematics.InverseKinematicsModel(
            robot=self.robot,
            iktype=IkParameterization.Type.Transform6D
            )
        
        if not  self.ikmodel.load():
             self.ikmodel.autogenerate()


        # Principal turbine elements linked to attributes
        self.bodies = self.env.GetBodies()
        
        try:
            self.rotor = next(body for body in self.bodies if body.GetName()==_ROTOR)
        except StopIteration:
            raise XMLStructError(_ROTOR)
        
        try:
            self.iris = next(body for body in self.bodies if body.GetName()==_IRIS)
        except StopIteration:
            raise XMLStructError(_IRIS)
        
        try:
            self.runner_area = next(body for body in self.bodies if body.GetName()==_RUNNER_AREA)
        except StopIteration:
            raise XMLStructError(_RUNNER_AREA)
        
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

            
    def place_rail(self,rail_place):
        """ Place both rails in the environment """
        #P - primary rail,
        #S - secondary rail,
        #alpha - angle from the perpendicular to the primary rail
        #Using camera standard axis and transformLookat for rails
        
        # Primary Rail
        primary_extent = self.primary.GetLinks()[0].GetGeometries()[0].GetBoxExtents()
        primary_offset = array([0,-primary_extent[1],-primary_extent[2]])+array([0, 0, self.config.environment.primary_safe_margin])
        primary_offset_transform = eye(4)
        primary_offset_transform[0:3,3] = primary_offset
        
        # Secondary Rail
        secondary_extent = self.secondary.GetLinks()[0].GetGeometries()[0].GetBoxExtents()
        #Resizing 
        secondary_extent[2] = abs(rail_place.s)/2.0 + self.config.environment.secondary_safe_margin
        self.env.RemoveKinBody(self.secondary)
        self.secondary.InitFromBoxes(array([concatenate([zeros(3),secondary_extent])]),True)
        self.env.AddKinBody(self.secondary)
        #
        secondary_offset = array([0,-secondary_extent[1],secondary_extent[2]])+array([0, self.config.environment.robot_level_difference, -self.config.environment.secondary_safe_margin])
        secondary_offset_transform = eye(4)
        secondary_offset_transform[0:3,3] = secondary_offset
        
        # Rails Traonsform and Placement
        primary_transform = transformLookat(array([0,0,self.config.environment.z_floor_level]),
                                            array([rail_place.p,0,self.config.environment.z_floor_level]),
                                            [0,0,1])
        self.primary.SetTransform(dot(primary_transform,primary_offset_transform))
        
        secondary_transform = transformLookat(rail_place.getXYZ(self),
                                                   array([rail_place.p,0,self.config.environment.z_floor_level]),
                                                   [0,0,1])    
        self.secondary.SetTransform(dot(secondary_transform,secondary_offset_transform))

    def check_rail_collision(self):
        """ Check rail <-> blades collision """
        with self.env:
            collisions = [self.env.CheckCollision(self.primary,blade)
                          or self.env.CheckCollision(self.secondary,blade)
                          for blade in self.blades]
        return any(collisions)

    def place_robot(self,rail_place):
        """ Place robot on the end of the secondary rail """
        Placement = eye(4)
        Placement[0:3,3] = rail_place.getXYZ(self) + [0, 0, self.config.environment.robot_level_difference]
        R = matrixFromAxisAngle([0, 0, rail_place.alpha])
        self.robot.SetTransform(dot(Placement, R))

    def check_robot_collision(self):
        """ Check robot <-> environment collision """
        with self.env:
            collisions = [self.env.CheckCollision(self.robot,body) for body in self.bodies]
        return any(collisions)
