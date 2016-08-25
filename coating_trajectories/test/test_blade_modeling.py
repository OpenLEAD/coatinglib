import os
from numpy import array, load, array_equal, abs, max, mean, sum, min

from . import TestCase
from .. import rbf
from ..turbine import Turbine

tolmax = 1e-2
tolmean = 1e-3
toleps = 1e-1

class TestBladeModeling(TestCase):
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
        cls.turb = Turbine(cls.test_dir + '/dummy.cfg',False)
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
        template_points = load(self.test_dir + '/template_points.npz')
        template_points = template_points['array']
        TestBladeModeling.blade1.sampling(delta = 0.005, min_distance_between_points=0.001)

        self.assertTrue(array_equal(TestBladeModeling.blade1._points,
                                    template_points))

    def test_make_model(self):
        template_points = load(self.test_dir + '/template_points.npz')
        template_points = template_points['array']
        TestBladeModeling.blade2._points = template_points

        TestBladeModeling.blade2.make_model()

        template_points = load(self.test_dir + '/template_extra_points.npz')
        template_points = template_points['array']

        rbf_results = []
        for point in template_points:
            rbf_results.append(TestBladeModeling.blade1._model.f(point[0:3]))
        rbf_results = abs(rbf_results)
        self.assertTrue(max(rbf_results)<tolmax and mean(rbf_results)<tolmean)

        outside_points = TestBladeModeling.blade2._model._pointsaugment(template_points)
        rbf_outside_results = []
        for point in outside_points:
            rbf_outside_results.append(TestBladeModeling.blade2._model.f(point[0:3]))
        rbf_outside_results = rbf_outside_results
        self.assertTrue(mean(rbf_outside_results)>TestBladeModeling.blade2._model._eps*(1-toleps),
                        msg = str(mean(rbf_outside_results)))

    def test_generate_trajectory(self):
        template_points = load(self.test_dir + '/template_points.npz')
        template_points = template_points['array']
        TestBladeModeling.blade3._points = template_points
        
        template_rbf_points = load(self.test_dir + '/template_rbf_points.npz')
        template_rbf_points = template_rbf_points['array']
        TestBladeModeling.blade3._model._points = template_rbf_points
        
        template_w = load(self.test_dir + '/template_w.npz')
        template_w = template_w['array']
        TestBladeModeling.blade3._model._w = template_w
        
        TestBladeModeling.blade3._modelLoaded = True

        template_trajectories = load(self.test_dir + '/template_trajectories.npz')
        template_trajectories = template_trajectories['array']

        sphere = mathtools.sphere(1.04, 0.97, 0.003)
        TestBladeModeling.blade3.generate_trajectory(sphere)

        for template_trajectory in template_trajectories:
            for point in template_trajectory:
                for trajectory in TestBladeModeling.blade3._trajectories:
                    dif = trajectory-point
                    if min(sum(dif*dif,1))<tolmean: break
                else: self.assertTrue(False)
        self.assertTrue(True)        
        
if __name__ == '__main__':
    unittest.main()
