from turbine import Turbine
from turbine_config import TurbineConfig, ConfigFileError
from visualizer import Visualizer
from os.path import join, realpath
from os import makedirs, environ
from openravepy import matrixFromAxisAngle, Sensor
from math import pi
import time
from copy import copy
from numpy import array

if __name__ == "__main__":
    
    dir_test = join(realpath('.'),'test')
    environ['OPENRAVE_DATA'] = str(dir_test)
    cfg = TurbineConfig.load('turbine_unittest.cfg','test')
    turb = Turbine(cfg)
    
    turb.env.Remove(turb.primary)
    turb.env.Remove(turb.secondary)
    turb.env.Remove(turb.iris)
    turb.env.Remove(turb.rotor)
    turb.env.Remove(turb.runner_area)
    turb.env.Remove(turb.blades[1])
    turb.env.Remove(turb.blades[2])
    turb.env.Remove(turb.blades[3])

    # Visualizer
    vis = Visualizer(turb.env)

    # Move Robot
    T = array([[-1,0,0,2.21],
               [0,0,1,-3.42],
               [0,1,0,0.69],
               [0,0,0,1]])
    turb.robot.SetTransform(T)

    env = turb.env
    sensor = env.GetSensors()[0]
    sensor.Configure(Sensor.ConfigureCommand.PowerOn)
    sensor.Configure(Sensor.ConfigureCommand.RenderDataOn)
    data = sensor.GetSensorData(Sensor.Type.Laser)
    dist = data.ranges[0]
    pos = data.positions[0]
    point = pos + dist
               
    
