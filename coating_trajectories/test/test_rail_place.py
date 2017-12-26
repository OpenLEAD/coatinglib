import unittest

from . import TestCase
from numpy import arange,argmax
from .. turbine_config import TurbineConfig, ConfigFileError
from .. import rail_place
from ..turbine import Turbine

class TestRailPlace(TestCase):
    @classmethod
    def setUpClass(cls):
        super(TestRailPlace, cls).setUpClass()
        turbconf = TurbineConfig.load("turbine_unittest.cfg", cls.test_dir)
        cls.turb = Turbine(turbconf)
        cls.samples =1000000
        cls.railplaces = rail_place.rand_rail(cls.turb.config, cls.samples)

    def test_railplace_constrains(self):
        total = 0
        alpha_min = TestRailPlace.turb.config.environment.rail_angle_mean - TestRailPlace.turb.config.environment.rail_angle_limit
        alpha_max = TestRailPlace.turb.config.environment.rail_angle_mean + TestRailPlace.turb.config.environment.rail_angle_limit
        for railpos in TestRailPlace.railplaces:
            x,y,z = railpos.getXYZ(TestRailPlace.turb.config)
            if ((TestRailPlace.turb.config.environment.x_min <= x <= TestRailPlace.turb.config.environment.x_max)
                and (TestRailPlace.turb.config.environment.y_min <= y <= TestRailPlace.turb.config.environment.y_max)
                and z==TestRailPlace.turb.config.environment.z_floor_level
                and (alpha_min <= railpos.alpha <= alpha_max)):
                total = total + 1
                
        self.assertEqual(total,TestRailPlace.samples, msg = "Apenas "+str(int(100.0*total/TestRailPlace.samples))+"% de acerto")

    def test_railplace_spread(self):
        total = 0
        alpha_min = TestRailPlace.turb.config.environment.rail_angle_mean - TestRailPlace.turb.config.environment.rail_angle_limit
        alpha_max = TestRailPlace.turb.config.environment.rail_angle_mean + TestRailPlace.turb.config.environment.rail_angle_limit
        railtuples = []
        for railpos in TestRailPlace.railplaces:
            x,y,_ = railpos.getXYZ(TestRailPlace.turb.config)
            railtuples.append((x,y,railpos.alpha))

        railtuples.sort()
        step = (TestRailPlace.turb.config.environment.x_max - TestRailPlace.turb.config.environment.x_min)*100.0/TestRailPlace.samples
        rng = arange(TestRailPlace.turb.config.environment.x_min+step,TestRailPlace.turb.config.environment.x_max+step,step)

        misses = 0
        total = len(rng)
        for xpos in rng:
            bins=[]
            for i in range(len(railtuples)):
                if railtuples[0][0] > xpos:
                    maxis = argmax(bins,axis=0)
                    if ((abs(bins[maxis[1]][1]-TestRailPlace.turb.config.environment.y_max)<=0.1*abs(TestRailPlace.turb.config.environment.y_min))
                        or (abs(bins[maxis[1]][2]-alpha_max)<=0.1*abs(alpha_max))
                        or (abs(bins[maxis[1]][2]-alpha_min)<=0.1*abs(alpha_min))):
                        misses = misses + 1
                    break
                bins.append(railtuples[0])
                del railtuples[0]

        self.assertAlmostEqual(1.0*misses/total,1.0,delta=0.05)
    
if __name__ == '__main__':
    unittest.main()
