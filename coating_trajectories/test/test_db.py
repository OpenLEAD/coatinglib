from numpy import array, dot, pi
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
        for k, v in ptn.iteritems():
            if v == 0:
                key = k
                break
        assert_allclose(key_00, key)

        key = 0
        ptn = self.DB_45.load_points_to_num()
        for k, v in ptn.iteritems():
            if v == 0:
                key = k
                break
        assert_allclose(key_45, key)


if __name__ == '__main__':
    unittest.main()
