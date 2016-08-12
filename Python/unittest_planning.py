import planning
import unittest
from numpy import linspace, cos, sin, arange, zeros, ones, random, array, array_equal, dstack
from math import pi
from copy import deepcopy

a = []
b = []

class TestPlanning(unittest.TestCase):

    def setUp(self):
        global a
        r = 1.43
        self.B = array([2.43, 0, -1])
        theta = linspace(-10,10,15)*pi/180
        x = cos(theta)*1.43
        y = sin(theta)*1.43
        xy = dstack((x,y))
        nxny = dstack((-y,x))
        z = arange(-0.05, 0.05, 0.003)
        nxnynz = dstack((nxny,zeros((len(x),))))
        N = len(theta)
        self.trajectories = []
        for zi in z:
            zr = ones((len(x),))*zi
            xyz = dstack((xy,zr))
            trajectory = zeros((N,6))
            trajectory[:,0:3] = xyz
            trajectory[:,3:6] = nxnynz
            self.trajectories.append(trajectory)
        a = deepcopy(self.trajectories)
        for traj in self.trajectories:
           random.shuffle(traj) 
             
    def test_sortTrajectories(self):
        global b
        b = planning.sortTrajectories2(self.trajectories)
        self.assertTrue(array_equal(self.orgtraj,
                                    planning.sortTrajectories2(self.trajectories)))             
            
if __name__ == '__main__':
    unittest.main()
