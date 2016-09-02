from numpy import array, load, array_equal, abs, max, mean, sum, min, vstack
from numpy import argmax, sign, zeros, dot, argmin, random, sqrt, concatenate

import unittest
from . import TestCase
from .. import rbf
from .. turbine import Turbine
from .. turbine_config import TurbineConfig, ConfigFileError
from .. import blade_modeling
from .. import mathtools

tolmax = 1e-2
tolmean = 1e-3
toleps = 1e-1

class TestBladeModeling(TestCase):
    
    @classmethod
    def setUpClass(cls):
        super(TestBladeModeling, cls).setUpClass()
        turbconf = TurbineConfig.load("/dummy.cfg", cls.test_dir)
        cls.turb = Turbine(turbconf, False)
        
    def setUp(self):
        name = "testblade"
        blade = TestBladeModeling.turb.blades[0]
        self.blade1 = blade_modeling.BladeModeling(name, TestBladeModeling.turb, blade)
        self.blade2 = blade_modeling.BladeModeling(name, TestBladeModeling.turb, blade)
        self.blade3 = blade_modeling.BladeModeling(name, TestBladeModeling.turb, blade)
        self.blade4 = blade_modeling.BladeModeling(name, TestBladeModeling.turb, blade)
        self.blade5 = blade_modeling.BladeModeling(name, TestBladeModeling.turb, blade)
        self.blade6 = blade_modeling.BladeModeling(name, TestBladeModeling.turb, blade)
        self.blade1 = blade_modeling.BladeModeling(name, TestBladeModeling.turb, blade)


    def test_sampling(self):
        """
        The sampling test uses the BladeModeling::sampling method to sample a cube with 100mm
        side length, which center is at (0,0,-1).
        """
        cube_edge = 0.05
        delta = 0.005
        self.blade1.sampling(delta, None)
        points_origin = self.blade1._points[:,0:3]

        # Checking if points are on the cube
        infinity_norm = abs(max(points_origin)-cube_edge)<=1e-4
        self.assertTrue(infinity_norm.all(),
                        msg="The point(s) is (are) not on the cube")

        # Checking if delta spacing between points is right
        delta_spacing = True
        for point in points_origin:
            dif = sum(abs(points_origin-point),1)
            dif = dif[dif>=1e-5] # excluding the evaluated point from list
            delta_spacing = delta_spacing and (abs(min(dif)-delta)<=1e-4)
            if ~delta_spacing:
                self.assertTrue(False, msg = 'point ='+str(point)+' is not delta spaced')

        # Checking if normal vectors are right
        points_origin_and_normals = zeros((len(points_origin),6))
        points_origin_and_normals[:,0:3] = points_origin
        points_origin_and_normals[:,3:6] = self.blade1._points[:,3:6]
        for point in points_origin_and_normals:
            abs_point = abs(point[0:3])
            if len(abs_point[abs(abs_point-max(abs_point))<=1e-5])>1:continue # Vertexes are exceptions
            else:
                normal = [0,0,0]
                normal[argmax(abs(point[0:3]))] = 1*sign(point[argmax(abs(point[0:3]))])
                cos_theta = dot(normal, point[3:6])
                if abs(cos_theta-1)>=1e-3:
                    self.assertTrue(False, msg = 'point ='+str(point)+' has wrong normal vector')           

    def test_filter_by_distance(self):
        """
        The test verifies if the distance between points are greater or equal a threshold.
        """
        threshold = 1
        points = random.uniform(-1,1,size=(100,6))
        self.blade6._points = points
        points = self.blade6.filter_by_distance(self.blade6._points, threshold)
        
        for point in points:
            dif = points[:,0:3]-point[0:3]
            euclidean_distance = sum(dif*dif,1)
            euclidean_distance = euclidean_distance[euclidean_distance>=1e-6] # excluding the evaluated point from list
            nearer_point = argmin(euclidean_distance)
            self.assertTrue(min(euclidean_distance)>=threshold**2,
                        msg = "The points: "+str(point)+" and "+str(points[nearer_point])+" are too close"
                            )

    def test_make_model(self):
        """
        The test verifies if the interpolation is right, using extra points of the cube (test object).
        Also, the it verifies outside points, and the normal vectors.
        """

        model = rbf.RBF('test','r3')
        delta = 0.008
        self.blade2.sampling(delta, None)
        template_points = self.blade2._points
        
        delta = 0.005
        self.blade2.sampling(delta, None)
        self.blade2.make_model(model)

        # Verifying if interpolation is right with the extra points of the cube 
        rbf_results = []
        for point in template_points:
            rbf_results.append(model.f(point[0:3]))
        rbf_results = abs(rbf_results)
        self.assertTrue(max(rbf_results)<tolmax and mean(rbf_results)<tolmean)

        # Verifying if the outside points have its mean approximately eps 
        outside_points = model._pointsaugment(template_points)
        rbf_outside_results = []
        for point in outside_points:
            rbf_outside_results.append(model.f(point[0:3]))
        rbf_outside_results = array(rbf_outside_results)
        self.assertTrue(abs(mean(rbf_outside_results)-model._eps)<=toleps,
                        msg = "The outside points mean equals "+str(mean(rbf_outside_results)))

        # Verifying if the normal vectors of extra points are right
        template_points[:,0:3] = template_points[:,0:3]
        for point in template_points:
            abs_point = abs(point[0:3])
            if len(abs_point[abs(abs_point-max(abs_point))<=1e-5])>1:continue # Vertexes are exceptions
            else:
                normal = [0,0,0]
                normal[argmax(abs(point[0:3]))] = 1*sign(point[argmax(abs(point[0:3]))])
                cos_theta = dot(normal, point[3:6])
                self.assertTrue(abs(cos_theta-1)<=1e-3, msg = 'point ='+str(point)+' has wrong normal vector')

    def test_compute_initial_point(self):
        """
        The test creates one hundred random sampling data with x^2+y^2+z^2 = 1, and random normal vectors.
        The point [sqrt(3),0,-1] with small disturbance is added to the data, and it will be the
        farest point, thus, it should be the initial point (answer of the test). 
        For the surfaces, the test creates a sphere with radius = 2, and center = (0,0,0), and a plane z=-1.
        The intersection between the surfaces is a known circle with radius = sqrt(3).
        The test verifies if the initial point is [sqrt(3),0,-1].
        Finally, the test verifies if, for a loaded trajectory with initial point [sqrt(3),0,-1], the
        next point is [sqrt(0.75), 0, -0.5] (the step of the sphere is 1 and the new plane is z=-0.5).
        """

        s = mathtools.Sphere(2, 0, 1)
        data = random.uniform(-1,1,size=(100,6))
        data[:,0:3] =  data[:,0:3]*(1.0/(sqrt(sum(data[:,0:3]*data[:,0:3],1)))).reshape(100,1)
        data[:,3:6] =  data[:,3:6]*(1.0/(sqrt(sum(data[:,3:6]*data[:,3:6],1)))).reshape(100,1)
        data[:,2]-=1
        right_initial_point = array([sqrt(3), 0, -1, 0, 0, 1])
        disturbance = array([float(random.rand(1)*1e-3), 0, 0, 0, 0, 0] )
        data = vstack((data, right_initial_point-disturbance))
        
        class ZPlane:
            def __init__(self, z0=-1):
                self.z = z0
            def f(self, p):
                return p[2]-self.z
            def df(self, p):
                return array([0, 0, 1])
        zp = ZPlane()
        
        self.blade3._points = data
        initial_point = self.blade3.compute_initial_point(zp, s)

        # Verify if the initial point is right
        dif = abs(initial_point - right_initial_point)<=1e-5
        self.assertTrue(dif.all(), msg = 'point ='+str(initial_point)+' is not the right initial point')

        # Verify if load trajectories and find next initial point is right
        self.blade3._trajectories = array([[right_initial_point]])
        zp = ZPlane(-0.5)
        initial_point = self.blade3.compute_initial_point(zp, s)
        right_initial_point = array([sqrt(0.75), 0, -0.5, 0, 0, 1])
        dif = abs(initial_point - right_initial_point)<=1e-5
        self.assertTrue(dif.all(), msg = 'point ='+str(initial_point)+' is not the right initial point')

    def test_draw_parallel(self):
        """
        The test creates a initial point [sqrt(3),0,-1]. For the surfaces, the test creates a sphere
        with radius = 2, and center = (0,0,0), and a plane z=-1. The intersection between the surfaces
        is the known circle with radius = sqrt(3). Therefore, all generated points should be on the
        implicit circle x^2+y^2 = 3, z = -1.
        """
        s = mathtools.Sphere(2, 0, 1)
        initial_point = array([sqrt(3), 0, -1, 0, 0, 1])
        
        class ZPlane:
            def __init__(self, z0=-1):
                self.z = z0
            def f(self, p):
                return p[2]-self.z
            def df(self, p):
                return array([0, 0, 1])
        zp = ZPlane()
        
        trajectory = self.blade4.draw_parallel(initial_point, zp, s, 0.001)
        trajectory = array(trajectory)

        # Verify if all points are on circle
        test_xy = abs(sum(trajectory[:,0:2]*trajectory[:,0:2],1)-3)<=1e-5 # verify x^2+y^2
        test_z = abs(trajectory[:,2]-zp.z)<=1e-5 # verify z
        test_circle = test_xy & test_z
        self.assertTrue(test_circle.all(), msg = 'point ='+str(trajectory[~test_circle])+' is not on trajectory')

    def test_generate_trajectories(self):
        """
        The test creates a initial point [sqrt(3),0,-1]. For the surfaces, the test creates an initial
        sphere with radius = 2, and center = (0,0,0), step = 0.5, stop = 0.9; and the model is a plane z=-z0.
        Therefore, the surfaces are:
        - sphere: radius = 2, model: z = -1, intersection: x^2+y^2 = 3
        - sphere: radius = 1.5, model: z = -0.75, intersection: x^2+y^2 = 1.6875
        - sphere: radius = 1.0, model: z = -0.5, intersection: x^2+y^2 = 0.75
        """
        s = mathtools.Sphere(2, 0.9, 0.5)
        initial_point = array([sqrt(3), 0, -1, 0, 0, 1])
        
        class ZPlane:
            def __init__(self, sphere, z0=-1):
                self.z = z0
                self._name = 'Plane'
                self.s = sphere
            def f(self, p):
                if abs(self.s._Rn-2)<=1e-5:
                    self.z = -1
                if abs(self.s._Rn-1.5)<=1e-5:
                    self.z = -0.75
                if abs(self.s._Rn-1)<=1e-5:
                    self.z = -0.5
                return p[2]-self.z
            def df(self, p):
                return array([0, 0, 1])
            def update(self,z):
                self.z = z
        zp = ZPlane(s)

        self.blade5._points = [initial_point]
        self.blade5.generate_trajectories(zp, s, 0.001)

        # Verify if all points are on circle 1
        trajectory = array(self.blade5._trajectories[0])
        test_xy = abs(sum(trajectory[:,0:2]*trajectory[:,0:2],1)-3)<=1e-5 # verify x^2+y^2
        test_z = abs(trajectory[:,2]+1)<=1e-5 # verify z
        test_circle = test_xy & test_z
        self.assertTrue(test_circle.all(), msg = 'point(s) ='+str(trajectory[~test_circle])+' is (are) not on trajectory')

        # Verify if all points are on circle 2
        trajectory = array(self.blade5._trajectories[1])
        test_xy = abs(sum(trajectory[:,0:2]*trajectory[:,0:2],1)-1.6875)<=1e-5 # verify x^2+y^2
        test_z = abs(trajectory[:,2]+0.75)<=1e-5 # verify z
        test_circle = test_xy & test_z
        self.assertTrue(test_circle.all(), msg = 'point(s) ='+str(trajectory[~test_circle])+' is (are) not on trajectory')

        # Verify if all points are on circle 3
        trajectory = array(self.blade5._trajectories[2])
        test_xy = abs(sum(trajectory[:,0:2]*trajectory[:,0:2],1)-0.75)<=1e-5 # verify x^2+y^2
        test_z = abs(trajectory[:,2]+0.5)<=1e-5 # verify z
        test_circle = test_xy & test_z
        self.assertTrue(test_circle.all(), msg = 'point(s) ='+str(trajectory[~test_circle])+' is (are) not on trajectory')
        
if __name__ == '__main__':
    unittest.main()
