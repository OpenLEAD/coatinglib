import unittest
import os

class TestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.test_dir = os.path.dirname(os.path.realpath(__file__))
        paths = ""
        if os.environ.has_key('OPENRAVE_DATA'):
            paths = os.environ['OPENRAVE_DATA']+":"
            
        os.environ['OPENRAVE_DATA'] = paths + str(cls.test_dir)
        super

    def setUp(self):
        self.test_dir = os.path.dirname(os.path.realpath(__file__))
        super
