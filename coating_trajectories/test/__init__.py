import unittest
import os

class TestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.test_dir = os.path.dirname(os.path.realpath(__file__))
        super

    def setUp(self):
        self.test_dir = os.path.dirname(os.path.realpath(__file__))
        super
