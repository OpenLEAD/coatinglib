import rbf
import blade_modeling
import math
import mathtools
from turbine import Turbine
import unittest
import os
from numpy import array, load, array_equal

tol = 1e-6

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
        
        rbf_model = rbf.RBF(name,'r3')
        cls.turb = Turbine('test/dummy.cfg',False)
        try:
            TestBladeModeling.turb.env.RemoveKinBody(TestBladeModeling.turb.primary)  
            TestBladeModeling.turb.env.RemoveKinBody(TestBladeModeling.turb.secondary)
            TestBladeModeling.turb.env.RemoveKinBody(TestBladeModeling.turb.robot)
            TestBladeModeling.turb.env.RemoveKinBody(TestBladeModeling.turb.runner_area)
            TestBladeModeling.turb.env.RemoveKinBody(TestBladeModeling.turb.iris)
        except: None
        cls.blade1 = blade_modeling.BladeModeling(name, rbf_model, TestBladeModeling.turb, False)
        cls.blade2 = blade_modeling.BladeModeling(name, rbf_model, TestBladeModeling.turb, False)
        cls.blade3 = blade_modeling.BladeModeling(name, rbf_model, TestBladeModeling.turb, False)


    def test_sampling(self):
        template_points = load('test/template_points.npz')
        template_points = template_points['array']
        TestBladeModeling.blade1.sampling(delta = 0.005, min_distance_between_points=0.001)
        self.assertTrue(array_equal(TestBladeModeling.blade1._points,
                                    template_points))

    def test_make_model(self):
        template_points = load('test/template_points.npz')
        template_points = template_points['array']
        TestBladeModeling.blade2._points = template_points

        template_rbf_points = load('test/template_rbf_points.npz')
        template_rbf_points = template_rbf_points['array']
        template_w = load('test/template_w.npz')
        template_w = template_w['array']

        
        TestBladeModeling.blade2.make_model()
        self.assertTrue(array_equal(TestBladeModeling.blade2._model._points,
                                    template_rbf_points))

        dif = abs(TestBladeModeling.blade2._model._w - template_w)
        
        self.assertTrue((dif<abs(tol*template_w)).all())

    def test_generate_trajectory(self):
        template_points = load('test/template_points.npz')
        template_points = template_points['array']
        TestBladeModeling.blade3._points = template_points
        
        template_rbf_points = load('test/template_rbf_points.npz')
        template_rbf_points = template_rbf_points['array']
        TestBladeModeling.blade3._model._points = template_rbf_points
        
        template_w = load('test/template_w.npz')
        template_w = template_w['array']
        TestBladeModeling.blade3._model._w = template_w
        
        TestBladeModeling.blade3._modelLoaded = True

        template_trajectories = load('test/template_trajectories.npz')
        template_trajectories = template_trajectories['array']

        sphere = mathtools.sphere(1.04, 0.97, 0.003)
        TestBladeModeling.blade3.generate_trajectory(sphere)

        dif = abs(TestBladeModeling.blade3._trajectories - template_trajectories)
        
        self.assertTrue((dif<abs(tol*template_trajectories)).all())
        

        

if __name__ == '__main__':
    unittest.main()
