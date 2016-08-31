import unittest
from numpy import array, array_equal, eye, pi, around, sqrt
from math import cos, sin

from .. import mathtools

class Testmathtools(unittest.TestCase):
    def test_hat(self):
        v = [1, 2, 3]
        M = array([[0, -3, 2],
                   [3, 0, -1],
                   [-2, 1, 0]])
        self.assertTrue(array_equal(M, mathtools.hat(v)))

    def test_Rab(self):
        R = array([[0, -1, 0],
                   [1, 0, 0],
                   [0, 0, 1]])
        self.assertTrue(array_equal(mathtools.Rab(array([1, 0, 0]), array([0, 1, 0])), R))

    def test_Raxis(self):
        M = array([[1, 0, 0],
                   [0, 0, -1],
                   [0, 1, 0]])
        self.assertTrue(array_equal(around(mathtools.Raxis(array([1, 0, 0]), pi/2),5), M))

    def test_surfaces_tangent(self):
        s = mathtools.sphere(Rn0=3)
        ray = array([3, 0, 0, 0, 0, 1])
        tan = mathtools.surfaces_tangent(ray, s)
        self.assertTrue(array_equal(array([0, 1, 0]), tan))

    def test_curvepoint(self):
        s = mathtools.Sphere(2, 0, 0.001)
        data = random.rand(100,6)
        for point in data:
            point[3:6] = point[3:6]/sqrt(dot(point[3:6],point[3:6]))
            point[0:3] = point[0:3]/sqrt(dot(point[0:3],point[0:3]))
        data[:,2]-=1
        
        class ZPlane:
            def __init__(self, z0=-1):
                self.z = z0
            def f(self, p):
                return p[2]-self.z
            def df(self, p):
                return array([0, 0, 1])
        zp = ZPlane()
        
        self.blade3._points = data
        self.blade3._model = zp
        initial_point = data[argmax(s.f_array(data))]
        s.findnextparallel(initial_point)
        initial_point = curvepoint(zp, s, p)

        # Testing if point belongs to the sphere and plane
        self.assertTrue(abs(s.f(initial_point))<=1e-5, msg = 'point ='+str(abs(s.f(initial_point)))+' does not belong to sphere')
        self.assertTrue(abs(zp.f(initial_point))<=1e-5, msg = 'point ='+str(abs(s.f(initial_point)))+' does not belong to plane')
        

        

if __name__ == '__main__':
    unittest.main()
