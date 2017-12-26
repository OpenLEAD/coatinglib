import unittest
from . import TestCase
from .. turbine import Turbine
from .. turbine_config import TurbineConfig
from .. import robot_utils
from numpy import concatenate, array, random, arange, dot
from numpy import sqrt, sin, cos, ones, sum, zeros, cross
from math import pi
from openravepy import RaveCreateCollisionChecker, CollisionOptions, CollisionReport, matrixFromAxisAngle
from .. import mathtools


class TestPlanning(TestCase):

    def setUp(self):
        self.path = "turbine_unittest.cfg"

    def test_ikfast(self):
        """
        The test verifies if the robot can reach a reachable point, a reachable point
        with the blade as an obstacle (flame can reach it as there is no collision),
        and a non reachable point. The RaveCreateCollisionChecker('ode') do not check
        collisions with an initial state, only timestep positions.
        """
        turbconf = TurbineConfig.load(self.path, self.test_dir)
        turb = Turbine(turbconf)
        turb.robot.GetLink('Flame').Enable(True)
        turb.env.SetCollisionChecker(RaveCreateCollisionChecker(turb.env,'ode'))


        # Reachable point
        point = turb.manipulator.GetTransform()[0:3,3]
        disturbance = random.uniform(-0.01,0.01,3)
        point = point + disturbance
        point = concatenate((point,-turb.manipulator.GetDirection()))

        iksol = robot_utils.ikfast(turb.robot, point)
        self.assertTrue(len(iksol)>0, msg = 'No solution was found')

        # Reachable point with blade as an obstacle
        blade_T = turb.robot.GetTransform()
        blade_T[0:3,3] = turb.manipulator.GetTransform()[0:3,3]+[0, 2, 0]
        turb.blades[0].SetTransform(blade_T)
        iksol = robot_utils.ikfast(turb.robot, point)
        self.assertTrue(len(iksol)>0, msg = 'No solution was found')

        # Non-reachable point
        point[0] += -10
        iksol = robot_utils.ikfast(turb.robot, point)
        self.assertTrue(len(iksol)==0, msg = 'A solution was found')

    def test_compute_robot_joints(self):
        """
        The test generates a trajectory and computes robot joints.
        """
        turbconf = TurbineConfig.load(self.path, self.test_dir)
        turb = Turbine(turbconf)
        robot_pos = turb.robot.GetTransform()[0:3,3]

        theta = arange(-pi/10,pi/10,0.001)
        # Generating positions
        x = -2.7 + (2 - cos(3*theta))*cos(theta)
        y = 1 + (2 - cos(3*theta))*sin(theta)
        z = zeros(len(theta))
        positions = [[x[i],y[i],z[i]] for i in range(0,len(x))]

        # Generating normal vectors
        ny = 2*cos(theta)+cos(2*theta)-2*cos(4*theta)
        nx = -2*sin(theta)+sin(2*theta)+2*sin(4*theta)
        nz = zeros(len(theta))
        normal_vectors = array([[nx[i],ny[i],nz[i]] for i in range(0,len(nx))])
        normal_vectors = normal_vectors*(1.0/sqrt(sum(normal_vectors*normal_vectors,1))).reshape(len(normal_vectors),1)
        for i in range(0,len(normal_vectors)):
            normal_vectors[i] = cross(normal_vectors[i],[0,0,1])

        trajectory = [[positions[i][0], positions[i][1], positions[i][2],
                       normal_vectors[i][0], normal_vectors[i][1],
                       normal_vectors[i][2]] for i in range(0,len(positions))]

        X90 = matrixFromAxisAngle([pi/2, 0, 0])
        Z180 = matrixFromAxisAngle([0, 0, pi])
        trajectory = array(mathtools.rotate_trajectories(turb, [trajectory], dot(Z180,X90))[0])
        trajectory[:,2] = trajectory[:,2] + robot_pos[2]

##        from .. import visualizer
##        vis = visualizer.Visualizer(turb.env)
##        vis.plot(array(trajectory)[:,0:3])
##        vis.plot_normal(trajectory)
##        x = raw_input('wait')

        joint_solutions = robot_utils.compute_robot_joints(turb, trajectory, 0)
        self.assertTrue(len(joint_solutions)==len(trajectory),
                        msg = 'The trajectory is not feasible')


if __name__ == '__main__':
    unittest.main()
