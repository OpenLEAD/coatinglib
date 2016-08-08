import RBF
import blade
import math
import mathtools
from turbine import Turbine
import unittest
import os

class TestBladeModeling(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(TestBladeModeling, cls).setUpClass()
        name = "testblade"
        
        try:
            os.remove('Blade/testblade_points.npz')
        except: None

        try:
            os.remove('Blade/Trajectory/testblade_r3_trajectories.npz')
        except: None

        try:
            os.remove('Blade/RBF/testblade_r3_points.npz')
        except: None

        try:
            os.remove('Blade/RBF/testblade_r3_w.npz')
        except: None
        
        rbf = RBF.RBF(name,'r3')
        cls.turb = Turbine('turbine_std.cfg',False)
        cls.blade = blade.BladeModeling(name, rbf, TestBladeModeling.turb, False)


    def test_modeling(self):
        self.assertTrue(TestBladeModeling.blade.sampling())
        self.assertTrue(TestBladeModeling.blade.make_model())
        sphere = mathtools.sphere(TestBladeModeling.turb.model.runner_radius,
                                  TestBladeModeling.turb.model.nose_radius,
                                  TestBladeModeling.turb.coating.parallel_gap
                                  )
        self.assertTrue(TestBladeModeling.blade.generate_trajectory(sphere))

        

if __name__ == '__main__':
    unittest.main()
