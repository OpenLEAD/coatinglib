import rbf
import blade_modeling
from turbine import Turbine
import mathtools
import unittest
import os
from numpy import array, load, abs, max, mean, sum, min

tolmax = 1e-2
tolmean = 1e-3
cosmean = 0.92

class TestRBF(unittest.TestCase):

    def setUp(self):

        self.template_points = load('test/template_extra_points.npz')
        self.template_points = self.template_points['array']
        self.template_normal = self.template_points[:,3:6]


    @classmethod
    def setUpClass(cls):
        super(TestRBF, cls).setUpClass()
        name = "test"
        template_points = load('test/template_points.npz')
        template_points = template_points['array']
     
        cls.rbf_r3 = rbf.RBF(name, 'r3', template_points)
        cls.rbf_logr = rbf.RBF(name, 'logr', template_points)
        cls.rbf_gaussr = rbf.RBF(name, 'gaussr', template_points)

    def test_r3(self):
        TestRBF.rbf_r3.make()

        rbf_results = []
        for point in self.template_points:
            rbf_results.append(TestRBF.rbf_r3.f(point[0:3]))

        rbf_results = abs(rbf_results)

        maxvalue = max(rbf_results)
        meanvalue = mean(rbf_results)

        self.assertTrue(maxvalue<tolmax and meanvalue<tolmean,
                        msg='max = '+str(maxvalue)+', mean = '+ str(meanvalue))

        rbf_normal = []
        for point in self.template_points:
            rbf_normal.append(TestRBF.rbf_r3.df(point[0:3]))
        rbf_normal = array(rbf_normal)
        cos_theta = sum(rbf_normal*self.template_normal,1)

        cos_theta_min = min(cos_theta)
        cos_theta_mean = mean(cos_theta)
        
        self.assertTrue(cos_theta_mean>cosmean,
                        msg='mean = '+ str(cos_theta_mean))
        

    def test_logr(self):
        TestRBF.rbf_logr.make()
       
        rbf_results = []
        for point in self.template_points:
            rbf_results.append(TestRBF.rbf_logr.f(point[0:3]))

        rbf_results = abs(rbf_results)

        maxvalue = max(rbf_results)
        meanvalue = mean(rbf_results)

        self.assertTrue(maxvalue<tolmax and meanvalue<tolmean,
                        msg='max = '+str(maxvalue)+', mean = '+ str(meanvalue))

        rbf_normal = []
        for point in self.template_points:
            rbf_normal.append(TestRBF.rbf_logr.df(point[0:3]))
        rbf_normal = array(rbf_normal)
        cos_theta = sum(rbf_normal*self.template_normal,1)

        cos_theta_min = min(cos_theta)
        cos_theta_mean = mean(cos_theta)
        
        self.assertTrue(cos_theta_mean>cosmean,
                        msg='mean = '+ str(cos_theta_mean))

    def test_gaussr(self):
        TestRBF.rbf_gaussr.make()

        rbf_results = []
        for point in self.template_points:
            rbf_results.append(TestRBF.rbf_gaussr.f(point[0:3]))

        rbf_results = abs(rbf_results)

        maxvalue = max(rbf_results)
        meanvalue = mean(rbf_results)

        self.assertTrue(maxvalue<tolmax and meanvalue<tolmean,
                        msg='max = '+str(maxvalue)+', mean = '+ str(meanvalue))

        rbf_normal = []
        for point in self.template_points:
            rbf_normal.append(TestRBF.rbf_gaussr.df(point[0:3]))
        rbf_normal = array(rbf_normal)
        cos_theta = sum(rbf_normal*self.template_normal,1)

        cos_theta_min = min(cos_theta)
        cos_theta_mean = mean(cos_theta)
        
        self.assertTrue(cos_theta_mean>cosmean,
                        msg='mean = '+ str(cos_theta_mean))

if __name__ == '__main__':
    unittest.main()
