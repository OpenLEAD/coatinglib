import os
from numpy import array, load, abs, max, mean, sum, min

import unittest
from . import TestCase
from .. import rbf

tolmax = 1e-2
tolmean = 1e-3
toleps = 1e-1
cosmean = 0.92

class TestRBF(TestCase):

    
    def setUp(self):
        template_points = load(self.test_dir + '/template/points.npz')
        template_points = template_points['array']
        
        name = "test"
        rbf_r3 = rbf.RBF(name, 'r3', template_points)
        rbf_logr = rbf.RBF(name, 'logr', template_points)
        rbf_gaussr = rbf.RBF(name, 'gaussr', template_points)
        
        self.template_points = load(self.test_dir + '/template/extra_points.npz')
        self.template_points = self.template_points['array']
        self.template_normal = self.template_points[:,3:6]

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

        outside_points = TestRBF.rbf_r3._pointsaugment(self.template_points)
        rbf_outside_results = []
        for point in outside_points:
            rbf_outside_results.append(TestRBF.rbf_r3.f(point[0:3]))
        rbf_outside_results = rbf_outside_results
        self.assertTrue(mean(rbf_outside_results)>TestRBF.rbf_r3._eps*(1-toleps),
                        msg = str(mean(rbf_outside_results)))

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

        outside_points = TestRBF.rbf_logr._pointsaugment(self.template_points)
        rbf_outside_results = []
        for point in outside_points:
            rbf_outside_results.append(TestRBF.rbf_logr.f(point[0:3]))
        rbf_outside_results = rbf_outside_results
        self.assertTrue(mean(rbf_outside_results)>TestRBF.rbf_logr._eps*(1-toleps),
                        msg = str(mean(rbf_outside_results)))

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

        outside_points = TestRBF.rbf_gaussr._pointsaugment(self.template_points)
        rbf_outside_results = []
        for point in outside_points:
            rbf_outside_results.append(TestRBF.rbf_gaussr.f(point[0:3]))
        rbf_outside_results = rbf_outside_results
        self.assertTrue(mean(rbf_outside_results)>TestRBF.rbf_gaussr._eps*(1-toleps),
                        msg = str(mean(rbf_outside_results)))

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
