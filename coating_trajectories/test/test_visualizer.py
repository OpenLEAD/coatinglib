from numpy import array, random, sum, sqrt

import unittest
from . import TestCase
from .. turbine import Turbine
from .. turbine_config import TurbineConfig, ConfigFileError
from .. visualizer import Visualizer

class TestVisualizer(TestCase):

    def setUp(self):
        turbconf = TurbineConfig.load("dummy.cfg", self.test_dir)
        turb = Turbine(turbconf)
        turb.env.Remove(turb.primary)
        turb.env.Remove(turb.secondary)
        turb.env.Remove(turb.iris)
        turb.env.Remove(turb.runner_area)
        turb.env.Remove(turb.robot)
        turb.env.Remove(turb.rotor)
        for i in range(0,len(turb.blades)):
            turb.env.Remove(turb.blades[i])
        
        self.vis = Visualizer(turb.env)

    def test_plot(self):
        """
        Verify if the visualizer can plot points, arrays, and arrays(arrays)
        """
        red = ((1,0,0))
        blue = ((0,0,1))
        green = ((0,1,0))
        
        # Test one point
        point = array([0,0,0])
        name = 'one_point'
        self.vis.plot(point, name, red, pointsize = 20)

        try:
            self.vis.handles['one_point']
            self.assertTrue(True)
        except KeyError:
            self.assertTrue(False, msg = 'point was not plotted')

        # Test array of points
        number_of_points = 10000
        points = random.uniform(-1,1,size=(number_of_points,3))
        points =  points*(1.0/(sqrt(sum(points*points,1)))).reshape(number_of_points,1)
        name = 'ten_thousand_points'
        self.vis.plot(points, name, blue)

        try:
            self.vis.handles[name]
            self.assertTrue(True)
        except KeyError:
            self.assertTrue(False, msg = 'points were not plotted')

        # Test array(array(points))
        number_of_points = 5000
        points_array = []
        number_of_arrays = 5
        for number in range(0,number_of_arrays):
            points = random.uniform(-1,1,size=(number_of_points,3))
            points =  points*(1.0/(sqrt(sum(points*points,1)))).reshape(number_of_points,1)
            center = random.uniform(-1,1,3)
            points = (points-center)*0.3
            points_array.append(points)
            
        name = 'one_thousand_array'
        self.vis.plot(points_array, name, green)

        try:
            self.vis.handles[name]
            self.assertTrue(True)
        except KeyError:
            self.assertTrue(False, msg = 'array of points was not plotted')

    
if __name__ == '__main__':
    unittest.main()
