import unittest
from . import TestCase
from .. turbine import Turbine
from .. turbine_config import TurbineConfig, ConfigFileError
from .. import planning
from numpy import concatenate

class TestPlanning(TestCase):

    def setUp(self):
        turbconf = TurbineConfig.load("/turbine_unittest.cfg", self.test_dir)
        turb = Turbine(turbconf, False)
        turb.env.Remove(turb.primary)
        turb.env.Remove(turb.secondary)
        turb.env.Remove(turb.iris)
        turb.env.Remove(turb.runner_area)
        turb.env.Remove(turb.robot)
        turb.env.Remove(turb.rotor)
        for i in range(0,len(turb.blades)):
            turb.env.Remove(turb.blades[i])
        self.turb = turb

    def test_ikfast(self):
        robot_transformation = self.turb.robot.GetTransform()
        point = robot_transformation[0:3,3]
        point = point + 0.5
        point = concatenate((point,[1,0,0]))
        facevector = [1,0,0]
        iksol = planning.ikfast(self.turb.ikmodel, facevector, point)
        print iksol


  
if __name__ == '__main__':
    unittest.main()
