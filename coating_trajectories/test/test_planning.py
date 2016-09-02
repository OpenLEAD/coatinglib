import unittest
from . import TestCase
from .. turbine import Turbine
from .. turbine_config import TurbineConfig, ConfigFileError
from .. import planning
from numpy import concatenate
from openravepy import IkParameterization

class TestPlanning(TestCase):

    def setUp(self):
        turbconf = TurbineConfig.load("/turbine_unittest.cfg", self.test_dir)
        turb = Turbine(turbconf, False)
        self.turb = turb

    def test_ikfast(self):
        self.turb.robot.GetLink('Flame').Enable(False)
        robot_transformation = self.turb.robot.GetTransform()
        point = robot_transformation[0:3,3]
        point = point + 0.5
        point = concatenate((point,[1,0,0]))
        facevector = [1,0,0]
        iksol = planning.ikfast(self.turb.ikmodel, facevector, point)
        self.assertTrue(len(iksol)>0, msg = 'no solution was found')


  
if __name__ == '__main__':
    unittest.main()
