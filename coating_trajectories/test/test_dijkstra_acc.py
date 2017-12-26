import unittest
from . import TestCase
from .. import dijkstra_acc
from numpy import array, load
from numpy import round
from os.path import join
import time


class TestPlanning(TestCase):

    def test_make_dijkstra(self):
        joints = [[[1.],[2.],[3.]],
                  [[10.],[20.],[30.],[40.]],
                  [[5.],[15.],[35.]],
                  [[100.]]]
        joints = [ array(j) for j in joints]

        dt = array([0.,5.,2.,10.])
        vel_limits = array([10.])
        acc_limits = array([1.])
        expected_path = [array([1.]),array([10.]),array([15.]),array([100.])]
        expected_min_cost = 0.8574
        joint_path, path, min_cost, adj, cost = dijkstra_acc.make_dijkstra(joints, dt, vel_limits, acc_limits, True)
        self.assertTrue(round(min_cost,5) == expected_min_cost, msg='min_cost is '+str(min_cost))

        for i in range(len(joint_path)):
            self.assertTrue(joint_path[i] == expected_path[i], msg='joint ' + str(joint_path[i]) + 'is not joint ' + str(expected_path[i]))

    def test_real_performance(self):
        joints = load(join('coating_trajectories','test','dijkstra_test','iksol.npz'))['array'][0]
        dt = load(join('coating_trajectories','test','dijkstra_test','dtimes.npz'))['array'][0]

        vel_limits = array([ 3.83972435, 3.4906585,  3.83972435, 7.15584993, 7.15584993])
        acc_limits = array([  6.70206433,  9.30260491, 15.98721595, 15.98721595, 15.98721595])
        t = time.time()
        dijkstra_acc.make_dijkstra(joints, dt, vel_limits, acc_limits, True)
        print '\nSingle parallel planing time is ', time.time() - t,'s. Target is less than 1s'


if __name__ == '__main__':
    unittest.main()