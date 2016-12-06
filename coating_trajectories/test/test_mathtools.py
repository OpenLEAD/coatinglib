import unittest
from numpy import array, array_equal, eye, pi, around, sqrt, random, dot
from math import cos, sin
from .. import mathtools
from . import TestCase

class Testmathtools(TestCase):

    def test_hat(self):
        v = [1, 2, 3]
        M = array([[0, -3, 2],
                   [3, 0, -1],
                   [-2, 1, 0]])
        self.assertTrue(array_equal(M, mathtools.hat(v)))

    def test_Rab(self):
        """
        Test the method for two random vectors, for known vectors and for the
        exception.
        """

        # Random vector test
        vector_1 = random.uniform(-1,1,3)
        vector_1 = vector_1/sqrt(dot(vector_1,vector_1))
        vector_2 = random.uniform(-1,1,3)
        vector_2 = vector_2/sqrt(dot(vector_2,vector_2))

        R = mathtools.Rab(vector_1, vector_2)
        v = dot(vector_2, R)
        self.assertTrue(max(abs(v-vector_1))<=1e-5)

        # Known vector test
        R = array([[0, -1, 0],
                   [1, 0, 0],
                   [0, 0, 1]])
        self.assertTrue(array_equal(mathtools.Rab(array([1, 0, 0]), array([0, 1, 0])), R))

        # Exception test
        vector_1 = [1,0,0]
        vector_2 = [-1,0,0]
        R = mathtools.Rab(vector_1, vector_2)
        v = dot(vector_2, R)
        self.assertTrue(max(abs(v-vector_1))<=1e-5)
        

    def test_compute_perpendicular_vector(self):
        """
        The test generates a random unit vector and verifies if a perpendicular
        vector is found with dot product.
        """
        vector_1  = random.uniform(-1,1,3)
        vector_1 = vector_1/sqrt(dot(vector_1,vector_1))
        perpendicular_vector = mathtools.compute_perpendicular_vector(vector_1)
        self.assertTrue(abs(dot(perpendicular_vector,vector_1))<=1e-6,
                        msg='vector_1:'+str(vector_1)+', perpendicular vector:'+str(perpendicular_vector))
    
    def test_Raxis(self):
        M = array([[1, 0, 0],
                   [0, 0, -1],
                   [0, 1, 0]])
        self.assertTrue(array_equal(around(mathtools.Raxis(array([1, 0, 0]), pi/2),5), M))

    def test_surfaces_tangent(self):
        s = mathtools.Sphere(Rn0=3)
        ray = array([3, 0, 0, 0, 0, 1])
        tan = mathtools.surfaces_tangent(ray, s)
        self.assertTrue(array_equal(array([0, 1, 0]), tan))

    def test_curvepoint(self):
        """
        The test creates an initial point with some disturbance.
        For the surfaces, the test creates a sphere with radius = 2, and center = (0,0,0); and a plane z=-1.
        The intersection between the surfaces is a known circle with radius = sqrt(3).
        The test verifies if the computed point belongs to sphere and plane.
        """
        s = mathtools.Sphere(2, 0, 0.001)
        initial_point = array([sqrt(3), 0, -1, 0, 0, 1])
        disturbance = array([float(random.rand(1)*1e-3), 0, 0, 0, 0, 0] )
        initial_point = initial_point + disturbance
        
        class ZPlane:
            def __init__(self, z0=-1):
                self.z = z0
            def f(self, p):
                return p[2]-self.z
            def df(self, p):
                return array([0, 0, 1])
        zp = ZPlane()
        
        s.findnextparallel(initial_point)
        initial_point = mathtools.curvepoint(zp, s, initial_point[0:3])

        # Testing if point belongs to the sphere and plane
        self.assertTrue(abs(s.f(initial_point))<=1e-5, msg = 'point ='+str(abs(s.f(initial_point)))+' does not belong to sphere')
        self.assertTrue(abs(zp.f(initial_point))<=1e-5, msg = 'point ='+str(abs(s.f(initial_point)))+' does not belong to plane')


    def test_filter_by_distance(self):
        """
        The test verifies if the distance between points are greater or equal a threshold.
        """

        threshold = 1
        points = random.uniform(-1,1,size=(100,6))
        points = mathtools.filter_by_distance(points, threshold)
        
        for point in points:
            dif = points[:,0:3]-point[0:3]
            euclidean_distance = sum(dif*dif,1)
            euclidean_distance = euclidean_distance[euclidean_distance>=1e-6] # excluding the evaluated point from list
            nearer_point = argmin(euclidean_distance)
            self.assertTrue(min(euclidean_distance)>=threshold**2,
                        msg = "The points: "+str(point)+" and "+str(points[nearer_point])+" are too close"
                            )



if __name__ == '__main__':
    unittest.main()
