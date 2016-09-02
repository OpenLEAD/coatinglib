import os
from numpy import array, load, abs, max, mean, sum, min, sqrt, random

import unittest
from . import TestCase
from .. import rbf

tolmax = 1e-2
tolmean = 1e-3
toleps = 1e-1
cosmean = 0.92

class TestRBF(TestCase):

    def setUp(self):
        """
        Making a sphere with random points for the RBFs.
        """
        number_of_model_data = 4000
        model_data = random.uniform(-1,1, size=(number_of_model_data,6))
        model_data[:,0:3] = model_data[:,0:3]*(1.0/(sqrt(sum(model_data[:,0:3]*model_data[:,0:3],1)))).reshape(number_of_model_data,1)
        model_data[:,3] = model_data[:,0]*2
        model_data[:,4] = model_data[:,1]*2
        model_data[:,5] = model_data[:,2]*2
        self.model_data = model_data

        number_of_validate_data = 100
        validate_data = random.uniform(-1,1, size=(number_of_validate_data,6))
        validate_data[:,0:3] = validate_data[:,0:3]*(1.0/(sqrt(sum(validate_data[:,0:3]*validate_data[:,0:3],1)))).reshape(number_of_validate_data,1)
        validate_data[:,3] = validate_data[:,0]*2
        validate_data[:,4] = validate_data[:,1]*2
        validate_data[:,5] = validate_data[:,2]*2
        self.validation_data = validate_data
        self.normal_validation_data = validate_data[:,3:6]

    def test_r3(self):
        """
        The test generates the RBF cubic kernel model with model_data samples,
        and computes the function for validation_data points. If mean and max values,
        outside points, and the normal vectors are in a threshold, the test passes. 
        """
        name = "test"
        rbf_r3 = rbf.RBF(name, 'r3', self.model_data)
        rbf_r3.make()

        # Verifying if validation_data is right
        rbf_results = []
        for point in self.validation_data:
            rbf_results.append(rbf_r3.f(point[0:3]))
        rbf_results = array(rbf_results)
        rbf_results = abs(rbf_results)
        maxvalue = max(rbf_results)
        meanvalue = mean(rbf_results)
        self.assertTrue(maxvalue<tolmax and meanvalue<tolmean,
                        msg = "max(validation_data) ="+str(maxvalue)+", it should be "+str(tolmax)+"mean(validation_data) ="+str(meanvalue)+", it should be "+str(tolmean))

        # Verifying if outside_points are right
        outside_points = rbf_r3._pointsaugment(self.validation_data)
        rbf_results = []
        for point in outside_points:
            rbf_results.append(rbf_r3.f(point[0:3]))
        rbf_results = array(rbf_results)
        self.assertTrue(abs(mean(rbf_results)-rbf_r3._eps)<=toleps,
                        msg = "mean(outside_points) ="+str(mean(rbf_results))+", it should be "+str(rbf_r3._eps))

        # Verifying if normal_vectors are right
        rbf_results = []
        for point in self.validation_data:
            rbf_results.append(rbf_r3.df(point[0:3]))
        rbf_results = array(rbf_results)
        cos_theta = sum(rbf_results*self.normal_validation_data,1)
        cos_theta_min = min(cos_theta)
        cos_theta_mean = mean(cos_theta)
        self.assertTrue(cos_theta_mean>cosmean,
                        msg='mean = '+ str(cos_theta_mean))
        
    def test_gaussr(self):
        """
        The test generates the RBF gaussian kernel model with model_data samples,
        and computes the function for validation_data points. If mean and max values,
        outside points, and the normal vectors are in a threshold, the test passes. 
        """
        name = "test"
        rbf_gaussr = rbf.RBF(name, 'gaussr', self.model_data)
        rbf_gaussr.make()

        # Verifying if validation_data is right
        rbf_results = []
        for point in self.validation_data:
            rbf_results.append(rbf_gaussr.f(point[0:3]))
        rbf_results = array(rbf_results)
        rbf_results = abs(rbf_results)
        maxvalue = max(rbf_results)
        meanvalue = mean(rbf_results)
        self.assertTrue(maxvalue<tolmax and meanvalue<tolmean,
                        msg = "max(validation_data) ="+str(maxvalue)+", it should be "+str(tolmax)+"mean(validation_data) ="+str(meanvalue)+", it should be "+str(tolmean))

        # Verifying if outside_points are right
        outside_points = rbf_gaussr._pointsaugment(self.validation_data)
        rbf_results = []
        for point in outside_points:
            rbf_results.append(rbf_gaussr.f(point[0:3]))
        rbf_results = array(rbf_results)
        self.assertTrue(abs(mean(rbf_results)-rbf_gaussr._eps)<=toleps,
                        msg = "mean(outside_points) ="+str(mean(rbf_results))+", it should be "+str(rbf_gaussr._eps))

        # Verifying if normal_vectors are right
        rbf_results = []
        for point in self.validation_data:
            rbf_results.append(rbf_gaussr.df(point[0:3]))
        rbf_results = array(rbf_results)
        cos_theta = sum(rbf_results*self.normal_validation_data,1)
        cos_theta_min = min(cos_theta)
        cos_theta_mean = mean(cos_theta)
        self.assertTrue(cos_theta_mean>cosmean,
                        msg='mean = '+ str(cos_theta_mean))

    def test_logr(self):
        """
        The test generates the RBF logarithmic kernel model with model_data samples,
        and computes the function for validation_data points. If mean and max values,
        outside points, and the normal vectors are in a threshold, the test passes. 
        """
        name = "test"
        rbf_logr = rbf.RBF(name, 'logr', self.model_data)
        rbf_logr.make()

        # Verifying if validation_data is right
        rbf_results = []
        for point in self.validation_data:
            rbf_results.append(rbf_logr.f(point[0:3]))
        rbf_results = array(rbf_results)
        rbf_results = abs(rbf_results)
        maxvalue = max(rbf_results)
        meanvalue = mean(rbf_results)
        self.assertTrue(maxvalue<tolmax and meanvalue<tolmean,
                        msg = "max(validation_data) ="+str(maxvalue)+", it should be "+str(tolmax)+"mean(validation_data) ="+str(meanvalue)+", it should be "+str(tolmean))

        # Verifying if outside_points are right
        outside_points = rbf_logr._pointsaugment(self.validation_data)
        rbf_results = []
        for point in outside_points:
            rbf_results.append(rbf_logr.f(point[0:3]))
        rbf_results = array(rbf_results)
        self.assertTrue(abs(mean(rbf_results)-rbf_logr._eps)<=toleps,
                        msg = "mean(outside_points) ="+str(mean(rbf_results))+", it should be "+str(rbf_logr._eps))

        # Verifying if normal_vectors are right
        rbf_results = []
        for point in self.validation_data:
            rbf_results.append(rbf_logr.df(point[0:3]))
        rbf_results = array(rbf_results)
        cos_theta = sum(rbf_results*self.normal_validation_data,1)
        cos_theta_min = min(cos_theta)
        cos_theta_mean = mean(cos_theta)
        self.assertTrue(cos_theta_mean>cosmean,
                        msg='mean = '+ str(cos_theta_mean))

if __name__ == '__main__':
    unittest.main()
