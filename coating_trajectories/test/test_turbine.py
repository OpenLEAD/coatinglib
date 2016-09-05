import unittest
from . import TestCase
from ..turbine_config import TurbineConfig, ConfigFileError
from ..turbine import Turbine

class TestTurbine(TestCase):

    def test_creation(self):
        turbconf = TurbineConfig.load("/turbine_unittest.cfg",TestTurbine.test_dir)
        turb = Turbine(turbconf)
        

if __name__ == '__main__':
    unittest.main()
