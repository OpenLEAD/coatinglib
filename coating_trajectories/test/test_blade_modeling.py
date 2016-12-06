from numpy import array, load, array_equal, abs, max, mean, sum, min, vstack
from numpy import argmax, sign, zeros, dot, argmin, random, sqrt, concatenate

import unittest
from . import TestCase
from .. import rbf
from .. turbine import Turbine
from .. turbine_config import TurbineConfig, ConfigFileError
from .. import blade_modeling
from .. import mathtools
import shutil
from copy import deepcopy, copy

tolmax = 1e-2
tolmean = 1e-3
toleps = 1e-1

class TestBladeModeling(TestCase):
    
    @classmethod
    def setUpClass(cls):
        super(TestBladeModeling, cls).setUpClass()
        turbconf = TurbineConfig.load("dummy.cfg", cls.test_dir)
        cls.turb = Turbine(turbconf)
        
    def setUp(self):
        self.name = "testblade"
        self.blade = TestBladeModeling.turb.blades[0]

    def test_sampling(self):
        """
        The sampling test uses the BladeModeling::sampling method to sample a cube with 100mm
        side length, which center is at (0,0,-1).
        """

        blade = blade_modeling.BladeModeling(TestBladeModeling.turb, self.blade)
        cube_edge = 0.05
        blade.samples_delta = 0.005
        delta = blade.samples_delta
        blade.min_distance_between_points = -1
        
        blade.sampling()
        points_origin = blade.points[:,0:3]

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
        points_origin_and_normals[:,3:6] = blade.points[:,3:6]
        for point in points_origin_and_normals:
            abs_point = abs(point[0:3])
            if len(abs_point[abs(abs_point-max(abs_point))<=1e-5])>1:continue # Vertexes are exceptions
            else:
                normal = [0,0,0]
                normal[argmax(abs(point[0:3]))] = 1*sign(point[argmax(abs(point[0:3]))])
                cos_theta = dot(normal, point[3:6])
                if abs(cos_theta-1)>=1e-3:
                    self.assertTrue(False, msg = 'point ='+str(point)+' has wrong normal vector')

        # Test saving samples and loading
        delta = copy(blade.samples_delta)
        shape = copy(blade.points.shape)
        min_distance_between_points = copy(blade.min_distance_between_points)
        blade.save_samples('test_sampling/', self.name)
        blade.load_samples('test_sampling/samples.xml')
        shutil.rmtree('test_sampling')
        
        self.assertTrue(abs(delta - blade.samples_delta)<=1e-5, msg = 'samples_delta was not loaded')
        self.assertTrue(abs(min_distance_between_points - blade.min_distance_between_points)<=1e-5,
                        msg = 'min_distance_between_points was not loaded')
        self.assertTrue(shape == blade.points.shape, msg = 'samples were not loaded')
        

    def test_make_model(self):
        """
        The test verifies if the interpolation is right, using extra points of the cube (test object).
        Also, it verifies outside points, and the normal vectors.
        """

        blade = blade_modeling.BladeModeling(TestBladeModeling.turb, self.blade)
        model = rbf.RBF('test','r3')
        blade.samples_delta = 0.008
        blade.min_distance_between_points = -1
        blade.sampling()
        template_points = blade.points

        blade.samples_delta = 0.005
        delta = blade.samples_delta
        blade.sampling()
        blade.make_model(model)
        model = blade.models[0]
        
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

        # Test save_model and loading
        model = deepcopy(blade.models[0])
        blade.save_model('test_make_model/', self.name)
        blade.load_model('test_make_model/model.xml')
        shutil.rmtree('test_make_model')

        self.assertTrue(blade.model_iter_surface is None, msg = 'model_iter_surface is not None, and it should be.')
        self.assertTrue(model.model_type == blade.models[0].model_type, msg = 'model_type check failed')
        for i in range(0,len(model._w)):
            self.assertTrue(abs(model._w[i]-blade.models[0]._w[i])<=1e-5, msg = 'model_w check failed')
            self.assertTrue((abs(model._points[i]-blade.models[0]._points[i])<=1e-5).all(), msg = 'model_point check failed')
        self.assertTrue(model._kernel == blade.models[0]._kernel, msg = 'kernel check failed')
        self.assertTrue(abs(model._eps-blade.models[0]._eps)<=1e-5, msg = 'eps check failed')
        

    def test_divide_model_points(self):
        """
        The test generates random data and call method for division.
        """

        blade = blade_modeling.BladeModeling(TestBladeModeling.turb, self.blade)
        number_of_data = 15000
        intersection_between_divisions = 0.15
        data = random.uniform(-10,10,size=(number_of_data,6))
        s = mathtools.Sphere(18, 0, 1)
        
        # Check output when number_of_points_per_part equals number of points
        model_points, model_index = blade.divide_model_points(data, s, number_of_data,
                                                              intersection_between_divisions)
        self.assertTrue(len(model_points)==1 and len(model_points[0])==len(data),
                        msg="Check fails when number_of_points_per_part equals number_of_data")

        # Check output when number_of_points_per_part is 1/3 of the number of points
        model_points, model_index = blade.divide_model_points(data, s, number_of_data/3,
                                                 intersection_between_divisions)

        self.assertTrue(len(model_points)==4,
                        msg="""Check fails when number_of_points_per_part is 1/3 of number_of_data
                            len(model_points) is not 4.    
                            """)
        check_model_points = True
        for i in range(0,len(model_points)):
            if i != len(model_points)-1:
                check_model_points = check_model_points and (len(model_points[i])==number_of_data/3)
        self.assertTrue(check_model_points,
                        msg="""Check fails when number_of_points_per_part is 1/3 of number_of_data
                            len(model_points[i]) is not number_of_data/3.    
                            """)

    def test_select_model(self):
        """
        The object is a cylinder. The zaxis is in (0,0). Data is generated as follows:
        Number of points per model = 2000 points
        Intersection coefficient = 0.15 (300 points)
        Groups: (0, 2.50), (2.00, 3.00), (2.70, 3.15), (3.10, 5.00)
        z = 0.0 .. 1.99 - 1700 points
        z = 2.00 .. 2.50 - 300 points
        z = 2.51 .. 2.69 - 1400 points
        z = 2.70 .. 3.00 - 300 points
        z = 3.01 .. 3.09 - 1400 points
        z = 3.10 .. 3.15 - 300 points
        z = 3.16 .. 5.00 - 1700 points
        R = 0.25
        A plane is the model_iter_surface.
        """

        blade = blade_modeling.BladeModeling(TestBladeModeling.turb, self.blade)
        model = rbf.RBF('test','r3')
        p = mathtools.Plane(0, 0, 1)
        blade.intersection_between_divisions = 0.15
        blade.number_of_points_per_model = 2000
        
        data = zeros((7100,6))
        xy = random.uniform(-1,1,size=(7100,2))
        R = (sqrt(sum(xy*xy,1))).reshape(7100,1)
        xy = xy/R
        xy = xy * sqrt(0.25)
        data[:,0:2] = xy[:,0:2]
        data[1700:2000, 2] = random.uniform(2, 2.5, 300)
        data[2000:3400, 2] = random.uniform(2.51, 2.69, 1400)
        data[3400:3700, 2] = random.uniform(2.70, 3.00, 300)
        data[3700:5100, 2] = random.uniform(3.01, 3.09, 1400)
        data[5100:5400, 2] = random.uniform(3.10, 3.15, 300)
        data[5400:7100, 2] = random.uniform(3.16, 5, 1700)
        data[:,3] = xy[:,0]
        data[:,4] = xy[:,1]

        blade.points = data
        blade.make_model(model, p)
        #blade.save_model('test_select_model/')

        test_data = zeros((8,6))
        xy = random.uniform(-1,1,size=(8,2))
        R = (sqrt(sum(xy*xy,1))).reshape(8,1)
        xy = xy/R
        xy = xy * sqrt(0.25)
        z = [1.5, 2.0, 2.51, 2.6,
             3.1, 3.12, 3.16, 4]
        test_data[:,0:2] = xy[:,0:2]
        test_data[:,3] = xy[:,0]
        test_data[:,4] = xy[:,1]
        test_data[:,2] = z

        self.assertTrue(len(blade.models)==4, msg='Number of models is'+str(len(blade.models)))
        index_verification = [0,0,1,1,2,2,3,3]

        for i in range(0,len(test_data)):
            model = blade.select_model(test_data[i])
            self.assertTrue(model == blade.models[index_verification[i]],
                            msg = 'Model was not well selected')
                

    def test_make_model_multiple(self):
        """
        The test verifies if the interpolation is right, using extra points
        and multiple RBF models.
        Tests the method select_model_from_list.
        """

        blade = blade_modeling.BladeModeling(TestBladeModeling.turb, self.blade)
        model = rbf.RBF('test','r3')
        s = mathtools.Sphere(0, 0, 1)
        
        number_of_model_data = 15000
        model_data = random.uniform(-1,1, size=(number_of_model_data,6))
        model_data[:,0:3] = model_data[:,0:3]*(1.0/(sqrt(sum(model_data[:,0:3]*model_data[:,0:3],1)))).reshape(number_of_model_data,1)
        model_data[:,2] = model_data[:,2]-1
        model_data[:,3] = model_data[:,0]*2
        model_data[:,4] = model_data[:,1]*2
        model_data[:,5] = model_data[:,2]*2

        number_of_validate_data = 100
        validate_data = random.uniform(-1,1, size=(number_of_validate_data,6))
        validate_data[:,0:3] = validate_data[:,0:3]*(1.0/(sqrt(sum(validate_data[:,0:3]*validate_data[:,0:3],1)))).reshape(number_of_validate_data,1)
        validate_data[:,2] = validate_data[:,2]-1
        validate_data[:,3] = validate_data[:,0]*2
        validate_data[:,4] = validate_data[:,1]*2
        validate_data[:,5] = validate_data[:,2]*2

        blade.points = model_data
        blade.number_of_points_per_model = number_of_model_data/3
        blade.intersection_between_divisions = 0.15
        
        blade.make_model(model, s)

        # Verifying if interpolation is right with extra points
        rbf_results = []
        for point in validate_data:
            model = blade.select_model(point)
            rbf_results.append(model.f(point))
        rbf_results = abs(rbf_results)
        self.assertTrue(max(rbf_results)<tolmax and mean(rbf_results)<tolmean)

        # Test save_model and loading
        models = deepcopy(blade.models)
        models_index = deepcopy(blade.models_index)
        model_iter_surface = deepcopy(blade.model_iter_surface)
        blade.save_model('test_make_model_multiple_RBF/', self.name)
        blade.load_model('test_make_model_multiple_RBF/model.xml')
        shutil.rmtree('test_make_model_multiple_RBF')

        self.assertTrue(model_iter_surface.name() == blade.model_iter_surface.name(), msg = 'model_iter_surface name check failed.')
        self.assertTrue(abs(model_iter_surface._Rn - blade.model_iter_surface._Rn)<=1e-5, msg = 'model_iter_surface Rn check failed')

        for i in range(0,len(blade.models)):
            for j in range(0,len(blade.models[i]._w)):
                self.assertTrue(abs(models[i]._w[j]-blade.models[i]._w[j])<=1e-5, msg = 'model_w check failed')
                self.assertTrue((abs(models[i]._points[j]-blade.models[i]._points[j])<=1e-5).all(), msg = 'model_point check failed')
            self.assertTrue(models[i]._kernel == blade.models[i]._kernel, msg = 'kernel check failed')
            self.assertTrue(abs(models[i]._eps-blade.models[i]._eps)<=1e-5, msg = 'eps check failed')
            self.assertTrue(abs(models_index[i][0]-blade.models_index[i][0])<=1e-5 and abs(models_index[i][1]-blade.models_index[i][1])<=1e-5,
            msg = 'models_index check failed')
            

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

        blade = blade_modeling.BladeModeling(TestBladeModeling.turb, self.blade)
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
        
        blade.points = data
        blade.models = [zp]
        initial_point = blade.compute_initial_point(s, [])

        # Verify if the initial point is right
        dif = abs(initial_point - right_initial_point)<=1e-5
        self.assertTrue(dif.all(), msg = 'point ='+str(initial_point)+' is not the right initial point')

        # Verify if load trajectories and find next initial point is right
        trajectories = array([[right_initial_point]])
        zp = ZPlane(-0.5)
        blade.models = [zp]
        initial_point = blade.compute_initial_point(s, trajectories)
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

        blade = blade_modeling.BladeModeling(TestBladeModeling.turb, self.blade)
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
        
        trajectory = blade.draw_parallel(initial_point, zp, s, 0.001)
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

        blade = blade_modeling.BladeModeling(TestBladeModeling.turb, self.blade)
        s = mathtools.Sphere(2, 0.9, 0.5)
        initial_point = array([sqrt(3), 0, -1, 0, 0, 1])
        
        class ZPlane:
            def __init__(self, sphere, z0=-1):
                self.z = z0
                self._name = 'Plane'
                self.s = sphere
                self.model_type = 'RBF'
                self._w = [0]
                self._eps = 0
                self._kernel = 'r3'
                self._points = [0]
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
        zp = ZPlane(s)

        blade.points = [initial_point]
        blade.models = [zp]
        blade.trajectory_step = 0.001
        trajectories = blade.generate_trajectories(s)

        # Verify if all points are on circle 1
        trajectory = array(trajectories[0])
        test_xy = abs(sum(trajectory[:,0:2]*trajectory[:,0:2],1)-3)<=1e-5 # verify x^2+y^2
        test_z = abs(trajectory[:,2]+1)<=1e-5 # verify z
        test_circle = test_xy & test_z
        self.assertTrue(test_circle.all(), msg = 'point(s) ='+str(trajectory[~test_circle])+' is (are) not on trajectory')

        # Verify if all points are on circle 2
        trajectory = array(trajectories[1])
        test_xy = abs(sum(trajectory[:,0:2]*trajectory[:,0:2],1)-1.6875)<=1e-5 # verify x^2+y^2
        test_z = abs(trajectory[:,2]+0.75)<=1e-5 # verify z
        test_circle = test_xy & test_z
        self.assertTrue(test_circle.all(), msg = 'point(s) ='+str(trajectory[~test_circle])+' is (are) not on trajectory')

        # Verify if all points are on circle 3
        trajectory = array(trajectories[2])
        test_xy = abs(sum(trajectory[:,0:2]*trajectory[:,0:2],1)-0.75)<=1e-5 # verify x^2+y^2
        test_z = abs(trajectory[:,2]+0.5)<=1e-5 # verify z
        test_circle = test_xy & test_z
        self.assertTrue(test_circle.all(), msg = 'point(s) ='+str(trajectory[~test_circle])+' is (are) not on trajectory')

        # Test save_trajectory
        trajectory_iter_surface = deepcopy(blade.trajectory_iter_surface)
        trajectory_step = blade.trajectory_step
        trajectories = deepcopy(blade.trajectories)
        
        blade.save_model('test/', self.name)
        blade.save_trajectory('test/model.xml', 'test_generate_trajectories/', self.name)
        blade.load_trajectory('test_generate_trajectories/trajectory.xml')
        shutil.rmtree('test_generate_trajectories')

        self.assertTrue(trajectory_iter_surface.name() == blade.trajectory_iter_surface.name(), msg = 'trajectory_iter_surface name check failed.')
        self.assertTrue(abs(trajectory_iter_surface._Rn0 - blade.trajectory_iter_surface._Rn)<=1e-5, msg = 'trajectory_iter_surface Rn check failed')
        self.assertTrue(abs(trajectory_iter_surface.stopR - blade.trajectory_iter_surface.stopR)<=1e-5, msg = 'trajectory_iter_surface stopR check failed')
        self.assertTrue(abs(trajectory_iter_surface.coatingstep - blade.trajectory_iter_surface.coatingstep)<=1e-5, msg = 'trajectory_iter_surface coatingstep check failed')

        for i in range(0, len(blade.trajectories)): # trajectory
            for j in range(0,len(trajectories[i])): # point
                self.assertTrue((abs(trajectories[i][j] - blade.trajectories[i][j])<=1e-5).all(), msg = 'trajectories check failed')
        
if __name__ == '__main__':
    unittest.main()
