import unittest
from . import TestCase
from .. turbine import Turbine
from .. turbine_config import TurbineConfig, ConfigFileError
from .. import dijkstra_vel
from numpy import concatenate, array, random, arange, dot, load
from numpy import sqrt, sin, cos, ones, sum, zeros, cross, round
from os.path import join
import time
from math import pi
from openravepy import RaveCreateCollisionChecker, CollisionOptions, CollisionReport, matrixFromAxisAngle
from .. import mathtools

class TestPlanning(TestCase):

    def setUp(self):
        self.path = "turbine_unittest.cfg"

#     def test_ikfast(self):
#         """
#         The test verifies if the robot can reach a reachable point, a reachable point
#         with the blade as an obstacle (flame can reach it as there is no collision),
#         and a non reachable point. The RaveCreateCollisionChecker('ode') do not check
#         collisions with an initial state, only timestep positions.
#         """
#         turbconf = TurbineConfig.load(self.path, self.test_dir)
#         turb = Turbine(turbconf)
#         turb.robot.GetLink('Flame').Enable(True)
#         turb.env.SetCollisionChecker(RaveCreateCollisionChecker(turb.env,'ode'))
#
#
#         # Reachable point
#         point = turb.manipulator.GetTransform()[0:3,3]
#         disturbance = random.uniform(-0.01,0.01,3)
#         point = point + disturbance
#         point = concatenate((point,-turb.manipulator.GetDirection()))
#
#         iksol = planning.ikfast(turb.robot, point)
#         self.assertTrue(len(iksol)>0, msg = 'No solution was found')
#
#         # Reachable point with blade as an obstacle
#         blade_T = turb.robot.GetTransform()
#         blade_T[0:3,3] = turb.manipulator.GetTransform()[0:3,3]+[0, 2, 0]
#         turb.blades[0].SetTransform(blade_T)
#         iksol = planning.ikfast(turb.robot, point)
#         self.assertTrue(len(iksol)>0, msg = 'No solution was found')
#
#         # Non-reachable point
#         point[0] += -10
#         iksol = planning.ikfast(turb.robot, point)
#         self.assertTrue(len(iksol)==0, msg = 'A solution was found')
#
#     def test_orientation_error_optimization(self):
#         """
#         Test a reachable point.
#         """
#         turbconf = TurbineConfig.load(self.path, self.test_dir)
#         turb = Turbine(turbconf)
#
#         lower_limits, upper_limits = turb.robot.GetActiveDOFLimits()
#         q = []
#         for joint_index in range(0,len(lower_limits)):
#             safe_region = (upper_limits[joint_index]+lower_limits[joint_index])/2
#             q.append(random.uniform(safe_region-safe_region*1e-5,
#                                     safe_region+safe_region*1e-5))
#
#         turb.robot.SetDOFValues(q)
#         Tee = turb.manipulator.GetTransform()
#         point_pos = Tee[0:3,3]
#         point_direction = -Tee[0:3,0]
#         point_direction = point_direction/sqrt(dot(point_direction,
#                                                    point_direction))
#         point = concatenate((point_pos, point_direction))
#         disturbance = random.uniform(-0.01,0.01,3)
#         point[0:3] = point[0:3]+disturbance
#         res = planning.orientation_error_optimization(turb, point)
#         self.assertTrue(res.success, msg = 'No solution was found')
#
#     def test_compute_robot_joints(self):
#         """
#         The test generates a trajectory and computes robot joints.
#         """
#         turbconf = TurbineConfig.load(self.path, self.test_dir)
#         turb = Turbine(turbconf)
#         robot_pos = turb.robot.GetTransform()[0:3,3]
#
#         theta = arange(-pi/10,pi/10,0.001)
#         # Generating positions
#         x = -2.7 + (2 - cos(3*theta))*cos(theta)
#         y = 1 + (2 - cos(3*theta))*sin(theta)
#         z = zeros(len(theta))
#         positions = [[x[i],y[i],z[i]] for i in range(0,len(x))]
#
#         # Generating normal vectors
#         ny = 2*cos(theta)+cos(2*theta)-2*cos(4*theta)
#         nx = -2*sin(theta)+sin(2*theta)+2*sin(4*theta)
#         nz = zeros(len(theta))
#         normal_vectors = array([[nx[i],ny[i],nz[i]] for i in range(0,len(nx))])
#         normal_vectors = normal_vectors*(1.0/sqrt(sum(normal_vectors*normal_vectors,1))).reshape(len(normal_vectors),1)
#         for i in range(0,len(normal_vectors)):
#             normal_vectors[i] = cross(normal_vectors[i],[0,0,1])
#
#         trajectory = [[positions[i][0], positions[i][1], positions[i][2],
#                        normal_vectors[i][0], normal_vectors[i][1],
#                        normal_vectors[i][2]] for i in range(0,len(positions))]
#
#         X90 = matrixFromAxisAngle([pi/2, 0, 0])
#         Z180 = matrixFromAxisAngle([0, 0, pi])
#         trajectory = array(mathtools.rotate_trajectories(turb, [trajectory], dot(Z180,X90))[0])
#         trajectory[:,2] = trajectory[:,2] + robot_pos[2]
#
# ##        from .. import visualizer
# ##        vis = visualizer.Visualizer(turb.env)
# ##        vis.plot(array(trajectory)[:,0:3])
# ##        vis.plot_normal(trajectory)
# ##        x = raw_input('wait')
#
#         joint_solutions = planning.compute_robot_joints(turb, trajectory, 0)
#         self.assertTrue(len(joint_solutions)==len(trajectory),
#                         msg = 'The trajectory is not feasible')


    def test_make_dijkstra(self):
        joints = [[[1.],[2.],[3.]],
                  [[10.],[20.],[30.],[40.]],
                  [[5.],[15.],[35.]],
                  [[100.]]]
        joints = [ array(j) for j in joints]

        dt = array([0.,5.,2.,10.])
        vel_limits = array([10.])
        acc_limits = array([1.])
        #expected_path = [array([1.]),array([10.]),array([15.]),array([100.])]
        expected_path = [array([3.]), array([10.]), array([15.]), array([100.])]
        expected_min_cost = 1.24
        #expected_min_cost = 0.8574
        joint_path, path, min_cost, adj, cost = dijkstra_vel.make_dijkstra_vel(joints, dt, vel_limits, True)
        self.assertTrue(round(min_cost,5) == expected_min_cost, msg='min_cost is '+str(min_cost))

        for i in range(len(joint_path)):
            self.assertTrue(joint_path[i] == expected_path[i], msg='joint ' + str(joint_path[i]) + 'is not joint ' + str(expected_path[i]))

    def test_real_dijkstra(self):
        joints = load(join('coating_trajectories','test','dijkstra_test','iksol.npz'))['array'][0]
        dt = load(join('coating_trajectories','test','dijkstra_test','dtimes.npz'))['array'][0]

        vel_limits = array([ 3.83972435, 3.4906585,  3.83972435, 7.15584993, 7.15584993])
        acc_limits = array([  6.70206433,  9.30260491, 15.98721595, 15.98721595, 15.98721595])
        expected_path = [array([1.]),array([10.]),array([15.]),array([100.])]
        expected_min_cost = 0.8574
        t = time.time()
        joint_path, path, min_cost, adj, cost = dijkstra_vel.make_dijkstra_vel(joints, dt, vel_limits, True)
        print 'Single parallel planing time is ', time.time() - t,'s. Target is less than 1s'
        # self.assertTrue(round(min_cost,5) == expected_min_cost, msg='min_cost is '+str(min_cost))
        #
        # for i in range(len(joint_path)):
        #     self.assertTrue(joint_path[i] == expected_path[i], msg='joint ' + str(joint_path[i]) + 'is not joint ' + str(expected_path[i]))

    def test_dijkstra_performance(self):
        n = 200
        joints = []
        for i in range(int(n)):
            joints.append(random.random((random.randint(1,30),1)))

        dt = random.random(len(joints))
        vel_limits = array([10.])
        acc_limits = array([1.])
        t = time.time()
        joint_path, path, min_cost, adj, cost = dijkstra_vel.make_dijkstra_vel(joints, dt, vel_limits, True)
        print 'time =', time.time() - t



if __name__ == '__main__':
    unittest.main()
