import unittest
from numpy import array, array_equal, pi, around, sqrt, random, dot, argmin, sin, cos, logspace, linalg, exp, linspace
from numpy.polynomial import legendre
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
        """ Test the method for two random vectors, for known vectors and for the
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
        """ The test generates a random unit vector and verifies if a perpendicular
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
        """ The test creates an initial point with some disturbance.
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
        """ The test verifies if the distance between points are greater or equal a threshold.
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

    def test_MLS(self):

        def asin(x):
            A = 2
            w = pi/4
            return A*sin(w*x)

        def dasin(x):
            A = 2
            w = pi / 4
            return A * w * cos(w * x)

        def ddasin(x):
            A = 2
            w = pi / 4
            return -A * w * w * sin(w * x)

        x = logspace(0,1,100)-1
        f = asin(x)
        df = dasin(x)
        ddf = ddasin(x)

        y, dy, ddy = mathtools.MLS(f,x,x,6)


        for i in range(len(y)):
            self.assertTrue(linalg.norm(f[i] - y[i]) <= 2e-3,
                            msg='position verification failed in ' + str(i) + '  ' + str(
                                linalg.norm(f[i] - y[i])) + ' ' + str(f[i]) + ',' + str(y[i]))
            self.assertTrue(linalg.norm(df[i] - dy[i]) <= 2e-3,
                            msg='joint velocity verification failed in ' + str(i) + '  ' + str(
                                linalg.norm(df[i] - dy[i])) + ' ' + str(df[i]) + ',' + str(dy[i]))
            self.assertTrue(linalg.norm(ddf[i] - ddy[i]) <= 1e-2,
                            msg='joint acc verification failed in ' + str(i) + '  ' + str(
                                linalg.norm(ddf[i] - ddy[i])) + ' ' + str(ddf[i]) + ',' + str(ddy[i]))

        def e(x):
            return exp(x)-1

        def de(x):
            return exp(x)

        def dde(x):
            return exp(x)

        x = linspace(0,4,100)
        f = e(x)
        df = de(x)
        ddf = dde(x)

        y, dy, ddy = mathtools.MLS(f,x,x,6)


        for i in range(len(y)):
            self.assertTrue(linalg.norm(f[i] - y[i]) <= 1e-2,
                            msg='position verification failed in ' + str(i) + '  ' + str(
                                linalg.norm(f[i] - y[i])) + ' ' + str(f[i]) + ',' + str(y[i]))
            self.assertTrue(linalg.norm(df[i] - dy[i]) <= 1e-1,
                            msg='joint velocity verification failed in ' + str(i) + '  ' + str(
                                linalg.norm(df[i] - dy[i])) + ' ' + str(df[i]) + ',' + str(dy[i]))
            self.assertTrue(linalg.norm(ddf[i] - ddy[i]) <= 1,
                            msg='joint acc verification failed in ' + str(i) + '  ' + str(
                                linalg.norm(ddf[i] - ddy[i])) + ' ' + str(ddf[i]) + ',' + str(ddy[i]))

    def test_legMLS(self):

        def asin(x):
            A = 2
            w = pi/4
            return A*sin(w*x)

        def dasin(x):
            A = 2
            w = pi / 4
            return A * w * cos(w * x)

        def ddasin(x):
            A = 2
            w = pi / 4
            return -A * w * w * sin(w * x)

        x = logspace(0,1,100)-1
        f = asin(x)
        df = dasin(x)
        ddf = ddasin(x)

        y, dy, ddy = mathtools.legMLS(f,x,x,6)


        for i in range(len(y)):
            self.assertTrue(linalg.norm(f[i] - y[i]) <= 1e-3,
                            msg='position verification failed in ' + str(i) + '  ' + str(
                                linalg.norm(f[i] - y[i])) + ' ' + str(f[i]) + ',' + str(y[i]))
            self.assertTrue(linalg.norm(df[i] - dy[i]) <= 1e-2,
                            msg='joint velocity verification failed in ' + str(i) + '  ' + str(
                                linalg.norm(df[i] - dy[i])) + ' ' + str(df[i]) + ',' + str(dy[i]))
            self.assertTrue(linalg.norm(ddf[i] - ddy[i]) <= 1e-1,
                            msg='joint acc verification failed in ' + str(i) + '  ' + str(
                                linalg.norm(ddf[i] - ddy[i])) + ' ' + str(ddf[i]) + ',' + str(ddy[i]))

        def e(x):
            return exp(x)-1

        def de(x):
            return exp(x)

        def dde(x):
            return exp(x)

        x = linspace(0,4,100)
        f = e(x)
        df = de(x)
        ddf = dde(x)

        y, dy, ddy = mathtools.legMLS(f,x,x,6)


        for i in range(len(y)):
            self.assertTrue(linalg.norm(f[i] - y[i]) <= 1e-3,
                            msg='position verification failed in ' + str(i) + '  ' + str(
                                linalg.norm(f[i] - y[i])) + ' ' + str(f[i]) + ',' + str(y[i]))
            self.assertTrue(linalg.norm(df[i] - dy[i]) <= 1e-2,
                            msg='joint velocity verification failed in ' + str(i) + '  ' + str(
                                linalg.norm(df[i] - dy[i])) + ' ' + str(df[i]) + ',' + str(dy[i]))
            self.assertTrue(linalg.norm(ddf[i] - ddy[i]) <= 1e-1,
                            msg='joint acc verification failed in ' + str(i) + '  ' + str(
                                linalg.norm(ddf[i] - ddy[i])) + ' ' + str(ddf[i]) + ',' + str(ddy[i]))

    def test_legcubic_path(self):

        def asin(x):
            A = 2
            w = pi/4
            return A*sin(w*x)

        def dasin(x):
            A = 2
            w = pi / 4
            return A * w * cos(w * x)

        def e(x):
            return exp(x)-1

        def de(x):
            return exp(x)

        p1 = array([asin(0), e(0)])
        p2 = array([asin(1), e(1)])
        dp1 = array([dasin(0), de(0)])
        dp2 = array([dasin(1), de(1)])

        c = mathtools.legcubic_path(p1,p2,dp1,dp2)
        legsol1 = legendre.legval(0, c)
        for i in range(len(p1)):
            self.assertAlmostEqual(legsol1[i],p1[i])

        legsol2 = legendre.legval(1, c)
        for i in range(len(p2)):
            self.assertAlmostEqual(legsol2[i], p2[i])

        legdsol1 = legendre.legval(0,legendre.legder(c))
        for i in range(len(dp1)):
            self.assertAlmostEqual(legdsol1[i], dp1[i])

        legdsol2 = legendre.legval(1, legendre.legder(c))
        for i in range(len(dp2)):
            self.assertAlmostEqual(legdsol2[i], dp2[i])

        self.assertEquals(legsol1.shape,p1.shape)
        self.assertEquals(legsol2.shape, p2.shape)
        self.assertEquals(legdsol1.shape, dp1.shape)
        self.assertEquals(legdsol2.shape, dp2.shape)

    def test_AccRampProfile(self):
        p0 = 5.
        dp0 = 1.
        dp1 = -2.5
        ddp0 = -3.
        ddp1 = 1.
        p1 = 7./24
        t1 = 1.
        t2= 2.

        accramp = mathtools.AccRampProfile(p0, p1, dp0, dp1, ddp0, ddp1, 2.5)

        self.assertAlmostEqual(ddp0, accramp.acc(0))
        self.assertAlmostEqual(ddp1, accramp.acc(2.5))
        self.assertAlmostEqual(dp0, accramp.vel(0))
        self.assertAlmostEqual(dp1, accramp.vel(2.5))
        self.assertAlmostEqual(p0, accramp.pos(0))
        self.assertAlmostEqual(p1, accramp.pos(2.5))

        self.assertAlmostEqual(t1, accramp.t1)
        self.assertAlmostEqual(t2, accramp.t2)

        self.assertAlmostEqual(-0.833333333333333, accramp.pos(3))
        self.assertAlmostEqual(4.5, accramp.pos(1))
        self.assertAlmostEqual(3.208333333333333, accramp.pos(1.5))
        self.assertAlmostEqual(-3.0, accramp.vel(1.5))

        accramp = mathtools.AccRampProfile(p0, p1, dp0, dp1, ddp0, ddp1, 0.)
        self.assertAlmostEqual(t1, accramp.t1)
        self.assertAlmostEqual(t2, accramp.t2)

    def test_AccDoubleStepProfile(self):
        p0 = 1; p1 = 2; v0 = -1; v1 = -2; a0 = -1; a1 = 2; amax = 7

        accdoublestep = mathtools.AccDoubleStepProfile(amax, p0, p1, v0, v1, a0, a1)








if __name__ == '__main__':
    unittest.main()
