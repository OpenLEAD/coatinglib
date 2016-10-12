import unittest
from . import TestCase
from .. turbine import Turbine
from .. turbine_config import TurbineConfig, ConfigFileError
from .. import planning
from numpy import concatenate, array, random, arange, dot
from numpy import sqrt, sin, cos, ones, sum, zeros, cross
from math import pi
from openravepy import RaveCreateCollisionChecker, CollisionOptions, CollisionReport

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
        point = concatenate((point,[1,0,0]))
        iksol = planning.ikfast(turb.robot, point)
        self.assertTrue(len(iksol)>0, msg = 'No solution was found')

        # Reachable point with blade as an obstacle
        blade_T = turb.robot.GetTransform()
        blade_T[0:3,3] = turb.manipulator.GetTransform()[0:3,3]+[0, 2, 0]
        turb.blades[0].SetTransform(blade_T)
        iksol = planning.ikfast(turb.robot, point)
        self.assertTrue(len(iksol)>0, msg = 'No solution was found')

        # Non-reachable point
        point[0] += -0.2
        iksol = planning.ikfast(turb.robot, point)
        self.assertTrue(len(iksol)==0, msg = 'A solution was found')

    def test_inverse_kinematics(self):
        """
        The test starts with a reachable point and modify it by a small step until the robot
        uses the tolerance to coat it. If there is no tolerance, the test will fail.
        """
        
        turbconf = TurbineConfig.load(self.path, self.test_dir)
        turb = Turbine(turbconf)
        
        step = 1e-3
        point = turb.manipulator.GetTransform()[0:3,3]
        point = concatenate((point,[1,0,0]))
        iksol,_ = planning.inverse_kinematics(turb, point)
        tolerance_verification = False

        while len(iksol)>0:
            turb.robot.SetDOFValues(iksol[0])
            point[0:3] = point[0:3]+step
            iksol, tolerance = planning.inverse_kinematics(turb, point)
            if tolerance:
                tolerance_verification = True
                break
        self.assertTrue(tolerance_verification, 'Tolerances could not be verified')

    def test_check_dof_limits(self):
        """
        Test a feasible and a non-feasible configuration of joint values. 
        """
        turbconf = TurbineConfig.load(self.path, self.test_dir)
        turb = Turbine(turbconf)
        robot = turb.robot
        q = robot.GetActiveDOFValues()
        self.assertTrue(planning.check_dof_limits(robot, q),
                        msg= 'DOF outside limits')

        q[0] = robot.GetJoints()[0].GetLimits()[0]
        self.assertTrue(~planning.check_dof_limits(robot, q),
                        msg='DOF inside limits')

    def test_orientation_error_optimization(self):
        """
        Test a reachable point. 
        """
        turbconf = TurbineConfig.load(self.path, self.test_dir)
        turb = Turbine(turbconf)

        lower_limits, upper_limits = turb.robot.GetActiveDOFLimits()
        q = []
        for joint_index in range(0,len(lower_limits)):
            safe_region = (upper_limits[joint_index]+lower_limits[joint_index])/2
            q.append(random.uniform(safe_region-safe_region*1e-5,
                                    safe_region+safe_region*1e-5))

        turb.robot.SetDOFValues(q)
        Tee = turb.manipulator.GetTransform()
        point_pos = Tee[0:3,3]
        point_direction = -Tee[0:3,0]
        point_direction = point_direction/sqrt(dot(point_direction,
                                                   point_direction))
        point = concatenate((point_pos, point_direction))
        disturbance = random.uniform(-0.01,0.01,3)
        point[0:3] = point[0:3]+disturbance
        res = planning.orientation_error_optimization(turb, point)
        self.assertTrue(res.success, msg = 'No solution was found')

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
        z = ones(len(theta))*robot_pos[2]
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

##        from .. import visualizer
##        vis = visualizer.Visualizer(turb.env)
##        vis.plot(array(trajectory)[:,0:3])
##        vis.plot_normal(trajectory)
##        x = raw_input('wait')
        joint_solutions = planning.compute_robot_joints(turb, trajectory, 0)
        self.assertTrue(len(joint_solutions)==len(trajectory),
                        msg = 'The trajectory is not feasible')
        


if __name__ == '__main__':
    unittest.main()
