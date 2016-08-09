import planning
import unittest
from numpy import linspace, cos, sin, arange, zeros, ones, random, array, array_equal
from math import pi
import copy

class TestPlanning(unittest.TestCase):

    def setUp(self):
        r = 1.43
        self.B = array([2.43, 0, -1])
        theta = linspace(-10,10,15)*pi/180
        x = cos(theta)*1.43
        y = sin(theta)*1.43
        z = arange(-0.05, 0.05, 0.003)
        nz = zeros((len(x),))
        N = len(theta)
        self.trajectories = []
        for zi in z:
            zr = ones((len(x),))*zi
            trajectory = zeros((N,6))
            trajectory[:,0] = x
            trajectory[:,1] = y
            trajectory[:,2] = zr
            trajectory[:,3] = -y
            trajectory[:,4] = x
            trajectory[:,5] = nz
            self.trajectories.append(trajectory)
        self.orgtraj = copy.copy(self.trajectories)
        for traj in self.trajectories:
           random.shuffle(traj) 
             
    def test_sortTrajectories(self):
        self.assertTrue(array_equal(self.trajectories,
                                    planning.sortTrajectories(self.B[0],
                                                              self.trajectories)))             
            
if __name__ == '__main__':
    unittest.main()
