from numpy import array, load, array_equal, abs, max, mean, sum, min

import unittest
from . import TestCase
from .. import rbf
from .. turbine import Turbine
from .. import blade_modeling
from .. import mathtools

tolmax = 1e-2
tolmean = 1e-3
toleps = 1e-1

class TestBladeModeling(TestCase):
    
    @classmethod
    def setUpClass(cls):
        super(TestBladeModeling, cls).setUpClass()
        cls.turb = Turbine("/dummy.cfg", cls.test_dir, False)
        
    def setUp(self):
        name = "testblade"
        rbf_model = rbf.RBF(name,'r3')
        self.blade1 = blade_modeling.BladeModeling(name, rbf_model, TestBladeModeling.turb)
        self.blade2 = blade_modeling.BladeModeling(name, rbf_model, TestBladeModeling.turb)
        self.blade3 = blade_modeling.BladeModeling(name, rbf_model, TestBladeModeling.turb)


    def test_sampling(self):
        template_points = load(self.test_dir + '/template/points.npz')
        template_points = template_points['array']
        self.blade1.sampling(delta = 0.005, min_distance_between_points=0.005)

        self.assertTrue(array_equal(self.blade1._points,
                                    template_points))

    def test_make_model(self):
        template_points = load(self.test_dir + '/template/points.npz')
        template_points = template_points['array']
        self.blade2._points = template_points

        self.blade2.make_model()

        template_points = load(self.test_dir + '/template/rbf_points.npz')
        template_points = template_points['array']

        rbf_results = []
        for point in template_points:
            rbf_results.append(self.blade1._model.f(point[0:3]))
        rbf_results = abs(rbf_results)
        self.assertTrue(max(rbf_results)<tolmax and mean(rbf_results)<tolmean)

        outside_points = self.blade2._model._pointsaugment(template_points)
        rbf_outside_results = []
        for point in outside_points:
            rbf_outside_results.append(self.blade2._model.f(point[0:3]))
        rbf_outside_results = rbf_outside_results
        self.assertTrue(mean(rbf_outside_results)>self.blade2._model._eps*(1-toleps),
                        msg = str(mean(rbf_outside_results)))

    def test_generate_trajectory(self):
        template_points = load(self.test_dir + '/template/points.npz')
        template_points = template_points['array']
        self.blade3._points = template_points
        
        template_rbf_points = load(self.test_dir + '/template/rbf_points.npz')
        template_rbf_points = template_rbf_points['array']
        self.blade3._model._points = template_rbf_points
        
        template_w = load(self.test_dir + '/template/w.npz')
        template_w = template_w['array']
        self.blade3._model._w = template_w

        template_trajectories = load(self.test_dir + '/template/trajectories.npz')
        template_trajectories = template_trajectories['array']

        step = 1e-3
        sphere = mathtools.Sphere(1.04, 0.97, 0.003)
        self.blade3.generate_trajectory(sphere, step)

        for template_trajectory in template_trajectories:
            for point in template_trajectory:
                for trajectory in self.blade3._trajectories:
                    dif = trajectory-point
                    if min(sum(dif*dif,1))<step*1e-1: break
                else: self.assertTrue(False)
        self.assertTrue(True)        
        
if __name__ == '__main__':
    unittest.main()
