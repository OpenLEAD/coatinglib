import mathtools
from turbine import Turbine
import unittest
from numpy import array, array_equal, eye, pi, around, sqrt
from math import cos, sin

class Testmathtools(unittest.TestCase):

    def setUp(self):
        self.turb = Turbine('turbine_std.cfg',False)
    
    def test_T(self):
        theta = pi/2

        Tx = array([[1, 0, 0, 0],
                   [0, cos(theta), -sin(theta), 0],
                   [0, sin(theta), cos(theta), 0],
                   [0, 0, 0, 1]])

        Ty = array([[cos(theta), 0, sin(theta), 0],
                   [0, 1, 0, 0],
                   [-sin(theta), 0, cos(theta), 0],
                   [0, 0, 0, 1]])

        Tz = array([[cos(theta), -sin(theta), 0, 0],
                   [sin(theta), cos(theta), 0, 0],
                   [0, 0, 1, 0],
                   [0, 0, 0, 1]])
        
        Te = eye(4)
        
        self.assertTrue(array_equal(mathtools.T(theta, 'x'), Tx))
        self.assertTrue(array_equal(mathtools.T(theta, 'y'), Ty))
        self.assertTrue(array_equal(mathtools.T(theta, 'z'), Tz))
        self.assertTrue(array_equal(mathtools.T(theta, ''), Te))

    def test_Rotate(self):
        obj = self.turb.blades[0]
        obj.SetTransform(eye(4))
        theta = pi/2

        Te = eye(4)

        obj = mathtools.Rotate(obj, -theta, 'z')
        obj = mathtools.Rotate(obj, theta, 'y')
        obj = mathtools.Rotate(obj, theta, 'z')
        obj = mathtools.Rotate(obj, theta, 'x')
        self.assertTrue(array_equal(obj.GetTransform(), Te))

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
        s1 = mathtools.sphere(Rn0=3)
        class plane:
            def __init__(self, y0=3.770):
                self.y = y0
            def f(self, p):
                return p[1]-self.y
            def df(self, p):
                return array([0, 1, 0])
        s2 = plane(0)    
        self.assertTrue(array_equal(array([0, 0, 3, 0, 0, 1]),
                                    around(mathtools.curvepoint(s1,
                                                                s2,
                                                                array([0, 0, 3.1])
                                                                )
                                           )
                                    )
                        )
        

if __name__ == '__main__':
    unittest.main()
