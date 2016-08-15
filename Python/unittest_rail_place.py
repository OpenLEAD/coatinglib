import unittest
import rail_place
from numpy import arange,argmax

class TestRailplace(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(TestRailplace, cls).setUpClass()
        cls.turb = rail_place.Turbine("turbine_unittest.cfg",False)
        cls.samples =1000000
        cls.railplaces = rail_place.rand_rail(cls.turb, cls.samples)

    def test_railplace_constrains(self):
        total = 0
        alpha_min = TestRailplace.turb.environment.rail_angle_mean - TestRailplace.turb.environment.rail_angle_limit
        alpha_max = TestRailplace.turb.environment.rail_angle_mean + TestRailplace.turb.environment.rail_angle_limit
        for railpos in TestRailplace.railplaces:
            x,y,z = railpos.getXYZ(TestRailplace.turb)
            if ((TestRailplace.turb.environment.x_min <= x <= TestRailplace.turb.environment.x_max)
                and (TestRailplace.turb.environment.y_min <= y <= TestRailplace.turb.environment.y_max)
                and z==TestRailplace.turb.environment.z_floor_level
                and (alpha_min <= railpos.alpha <= alpha_max)):
                total = total + 1
                
        self.assertEqual(total,TestRailplace.samples, msg = "Apenas "+str(int(100.0*total/TestRailplace.samples))+"% de acerto")

    def test_railplace_spread(self):
        total = 0
        alpha_min = TestRailplace.turb.environment.rail_angle_mean - TestRailplace.turb.environment.rail_angle_limit
        alpha_max = TestRailplace.turb.environment.rail_angle_mean + TestRailplace.turb.environment.rail_angle_limit
        railtuples = []
        for railpos in TestRailplace.railplaces:
            x,y,_ = railpos.getXYZ(TestRailplace.turb)
            railtuples.append((x,y,railpos.alpha))

        railtuples.sort()
        step = (TestRailplace.turb.environment.x_max - TestRailplace.turb.environment.x_min)*100.0/TestRailplace.samples
        rng = arange(TestRailplace.turb.environment.x_min+step,TestRailplace.turb.environment.x_max+step,step)

        misses = 0
        total = len(rng)
        for xpos in rng:
            bins=[]
            for i in range(len(railtuples)):
                if railtuples[0][0] > xpos:
                    maxis = argmax(bins,axis=0)
                    if ((abs(bins[maxis[1]][1]-TestRailplace.turb.environment.y_max)<=0.1*abs(TestRailplace.turb.environment.y_min))
                        or (abs(bins[maxis[1]][2]-alpha_max)<=0.1*abs(alpha_max))
                        or (abs(bins[maxis[1]][2]-alpha_min)<=0.1*abs(alpha_min))):
                        misses = misses + 1
                    break
                bins.append(railtuples[0])
                del railtuples[0]

        self.assertAlmostEqual(1.0*misses/total,1.0,delta=0.05)
    
if __name__ == '__main__':
    unittest.main()
