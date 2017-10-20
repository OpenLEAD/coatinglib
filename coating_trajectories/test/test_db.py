from numpy import array, dot, pi, concatenate
from numpy.testing import assert_allclose
from openravepy import matrixFromAxisAngle
import unittest
from . import TestCase
from .. import db
from ..turbine_config import TurbineConfig
from ..turbine import Turbine

class Testdb(TestCase):
    def setUp(self):
        turbconf = TurbineConfig.load("turbine_unittest.cfg", "")
        turbine = Turbine(turbconf)
        self.DB_00 = db.DB('DB_TEST', turbine, 'db')
        self.DB_45 = db.DB('DB_TEST', turbine, 'db_45')

        try:
            self.DB_m45 = db.DB('DB_TEST', turbine, 'db_-45')
            raise NameError("Unexpected pass")
        except:
            None

        self.DB = db.DB('DB_TEST', turbine, '')
        self.DB_VERT = db.DB('DB_TEST2', turbine, 'db_00')

    def test_init(self):
        self.assertTrue(self.DB_00.db_main_path == 'DB_TEST/db_00', msg='wrong path')
        self.assertTrue(self.DB_45.db_main_path == 'DB_TEST/db_45', msg='wrong path')
        self.assertTrue(self.DB.db_main_path == 'DB_TEST', msg='wrong path')
        self.assertTrue(self.DB_VERT.db_main_path == 'DB_TEST2/db_00', msg='wrong path')

        self.assertTrue(self.DB_00.verticalized == 0, msg='db is not verticalized')
        self.assertTrue(self.DB_45.verticalized == 0, msg='db is not verticalized')
        self.assertTrue(self.DB.verticalized == 0, msg='db is not verticalized')
        self.assertTrue(self.DB_VERT.verticalized == 1, msg='db is verticalidez')

        T = array([[1., 0., 0., 0.],
                   [0., 1., 0., 0.],
                   [0., 0., 1., 0.],
                   [0., 0., 0., 1.]])
        T_00 = array([[1., 0., 0., 0.],
                      [0., 0., 1., 0.],
                      [0., -1., 0., 0.],
                      [0., 0., 0., 1.]])
        T_VERT = T_00
        T_45 = array([[1., 0., 0., 0.],
                      [0., 0.70710678, 0.70710678, 0.],
                      [0., -0.70710678, 0.70710678, 0.],
                      [0., 0., 0., 1.]])

        assert_allclose(T, self.DB.T)
        assert_allclose(T_00, self.DB_00.T)
        assert_allclose(T_VERT, self.DB_VERT.T)
        assert_allclose(T_45, self.DB_45.T)

    def test_load_db(self):
        main_db = self.DB_00.load()
        assert_allclose(array(list(main_db[215])), array([205, 87]))

    def test_load_points_to_num(self):
        key_00 = array((-0.152195798, 0.538141088, -3.765726656))
        T = matrixFromAxisAngle([pi/4,0,0])
        key_45 = dot(T,key_00)

        key = 0
        ptn = self.DB_00.load_points_to_num()
        self.assertTrue(isinstance(ptn, dict), msg='points_to_num of DB_00 should be a dictionary')
        for k, v in ptn.iteritems():
            if v == 0:
                key = k
                break
        assert_allclose(key_00, key)

        key = 0
        ptn_45 = self.DB_45.load_points_to_num()
        self.assertTrue(isinstance(ptn_45, dict), msg='points_to_num of DB_45 should be a dictionary')
        for k, v in ptn_45.iteritems():
            if v == 0:
                key = k
                break
        assert_allclose(key_45, key)

    def test_grid_to_trajectories(self):

        points = array([2831, 2832, 2833, 2777, 2778, 2779, 2780])
        borders = array([[-0.34395401, 0.78002898, -3.65800373], [0.04986939, -0.00596685, -3.75569145]])
        point_vert = array([-0.94672199,  0.76057971, -1.52572671, -0.82705471, -0.51077813,
       -0.23470453])
        T = matrixFromAxisAngle([pi / 4, 0, 0])

        gtt_00 = self.DB_00.load_grid_to_trajectories()
        self.assertTrue(isinstance(gtt_00,dict), msg='grid_to_trajectories 00 should be a dictionary')
        self.assertTrue(len(gtt_00[0])==2,
                        msg='grid_to_trajectories[0] should be a size 2 list (trajectories and borders)')
        self.assertTrue(len(gtt_00[0][0]) == 151,
                        msg='grid_to_trajectories[0][0] should have 151 trajectories')
        assert_allclose(gtt_00[0][0][0], points)
        assert_allclose(gtt_00[0][1][0], borders)
        borders_45 = [dot(T[0:3,0:3],gtt_00[0][1][0][0]),dot(T[0:3,0:3],gtt_00[0][1][0][1])]

        gtt_45 = self.DB_45.load_grid_to_trajectories()
        self.assertTrue(isinstance(gtt_45, dict), msg='grid_to_trajectories 45 should be a dictionary')
        self.assertTrue(len(gtt_45[0]) == 2,
                        msg='grid_to_trajectories[0] should be a size 2 list (trajectories and borders)')
        self.assertTrue(len(gtt_45[0][0]) == 151,
                        msg='grid_to_trajectories[0][0] should have 151 trajectories')
        assert_allclose(gtt_45[0][0][0], points)
        assert_allclose(gtt_45[0][1][0], borders_45)

        gtt_vert = self.DB_VERT.load_grid_to_trajectories()
        self.assertTrue(isinstance(gtt_vert, dict), msg='grid_to_trajectories vert should be a dictionary')
        self.assertTrue(len(gtt_vert[64]) == 26,
                        msg='grid_to_trajectories[0] should be a size 26 list (trajectories)')
        self.assertTrue(len(gtt_vert[64][0]) == 118,
                        msg='grid_to_trajectories[0][0] should have 118 points')
        assert_allclose(gtt_vert[64][0][0], point_vert)

    def test_load_blade(self):
        point_00 = array([-0.1521958 ,  0.53814109, -3.76572666,  0.08927295,  0.16974092,
       -0.98143689])
        T = matrixFromAxisAngle([pi / 4, 0, 0])

        blade_00 = self.DB.load_blade()
        self.assertTrue(isinstance(blade_00.trajectories,list))
        assert_allclose(blade_00.trajectories[0][0],point_00)
        point_45 = concatenate((dot(T[0:3, 0:3], blade_00.trajectories[0][0][0:3]),
                                dot(T[0:3, 0:3], blade_00.trajectories[0][0][3:6])))

        blade_45 = self.DB_45.load_blade()
        self.assertTrue(isinstance(blade_45.trajectories, list))
        assert_allclose(blade_45.trajectories[0][0], point_45)

    def test_get_sorted_bases(self):
        ntb_0 = (1.1119167721707579, -0.45288596503722783, 0.24255236565450133)
        ntb_00 = self.DB_00.get_sorted_bases()
        ntb_45 = self.DB_45.get_sorted_bases()
        assert_allclose(ntb_00[0], ntb_0)
        assert_allclose(ntb_45[0], ntb_0)

    def test_get_sorted_points(self):
        point_0 = (-0.152195798, 0.538141088, -3.765726656)
        T = matrixFromAxisAngle([pi / 4, 0, 0])
        ntp_00 = self.DB_00.get_sorted_points()
        assert_allclose(ntp_00[0], point_0)
        point_45 = dot(T[0:3,0:3], point_0)

        ntp_45 = self.DB_45.get_sorted_points()
        assert_allclose(ntp_45[0], point_45)




if __name__ == '__main__':
    unittest.main()
