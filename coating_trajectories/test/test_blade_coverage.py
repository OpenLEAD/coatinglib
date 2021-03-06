import unittest
from . import TestCase
from .. turbine import Turbine
from .. turbine_config import TurbineConfig
from numpy import zeros, linspace, linalg, array, dot, logspace, hstack, einsum, ones, cumsum, vstack, abs
from math import cos, pi
from .. import blade_coverage
from openravepy import ConfigurationSpecification, RaveCreateTrajectory
import shutil
from .. import mathtools
from scipy.integrate import ode

## @file
# @brief This is a unittest for blade_coverage
# @author Renan S. Freitas
# @bug No known bugs

class TestBladeCoverage(TestCase):
    @classmethod
    def setUpClass(cls):
        super(TestBladeCoverage, cls).setUpClass()
        turbconf = TurbineConfig.load("dummy.cfg", cls.test_dir)
        cls.turbine = Turbine(turbconf)

    def setUp(self):
        """ This setUp creates reachable waypoints (rays) by moving the robot from (minimum joint limits)/10 to
        (maximum joint limits)/10 and storing the end-effector positions-orientations.
        """

        self.robot = TestBladeCoverage.turbine.robot
        minL, maxL = self.robot.GetDOFLimits()
        self.rays = zeros((100,6))
        self.joints = zeros((100,6))
        for i in range(len(minL)):
            self.joints[:, i] = linspace(minL[i] / 10.0, maxL[i] / 10.0, 100)
        with self.robot:
            for i, joint in enumerate(self.joints):
                self.robot.SetDOFValues(joint)
                T = self.robot.GetActiveManipulator().GetTransform()
                self.rays[i][0:3] = T[0:3,3]
                n = T[0:3,0]
                n = n/linalg.norm(n)
                self.rays[i][3:6] = -n
        self.rays = [self.rays, array(list(reversed(self.rays)))]
        self.turbine = TestBladeCoverage.turbine
        self.threshold = 5e-2

    def test_execute(self):
        """ This test verifies if the robot can complete the path. The waypoints were generated in the setUp and are
         feasible.
        """

        path = blade_coverage.Path(self.rays)
        path.execute(self.turbine, self.threshold)
        self.assertTrue(path.success, msg='path was not successfully executed')

    def test_serialize(self):
        """ This test verifies if the class can serialize a path in openrave format.
        """

        ind = str()
        for i in range(self.robot.GetDOF()): ind += str(i) + ' '

        cs = ConfigurationSpecification()
        _ = cs.AddGroup('joint_values ' + self.robot.GetName() + ' ' + ind, len(self.robot.GetActiveDOFIndices()),
                        'cubic')

        traj = RaveCreateTrajectory(self.turbine.env, '')
        traj.Init(cs)
        for i in range(len(self.joints)):
            waypoint = list(self.joints[i])
            traj.Insert(traj.GetNumWaypoints(), waypoint)
        path = blade_coverage.Path(self.rays)
        path.data = [traj]
        try:
            path.serialize('test_serialize')
            shutil.rmtree('test_serialize')
        except:
            self.assertTrue(False, msg='Data can not be serialized')

    def test_deserialize(self):
        """ This test verifies if the class can deserialize a file path in openrave format.
        """

        path = blade_coverage.Path(self.rays)
        try:
            path.deserialize(self.turbine, 'coating_trajectories/test/test_deserialize')
        except:
            self.assertTrue(False, msg='Data can not be deserialized')

    def test_get_joint(self):
        """ This test creates a path in openrave format and test the class method get_joint to return a known joint.
        """

        path = blade_coverage.Path(self.rays)
        ind = str()
        for i in range(self.robot.GetDOF()): ind += str(i) + ' '

        cs = ConfigurationSpecification()
        _ = cs.AddGroup('joint_values ' + self.robot.GetName() + ' ' + ind, len(self.robot.GetActiveDOFIndices()),
                        'cubic')
        cs.AddDerivativeGroups(1, False)
        cs.AddDerivativeGroups(2, False)
        _ = cs.AddDeltaTimeGroup()

        traj = RaveCreateTrajectory(self.turbine.env, '')
        traj.Init(cs)
        for i in range(len(self.joints[:1])):
            waypoint = list(self.joints[i])
            waypoint.extend(list(self.joints[i]))
            waypoint.extend(list(self.joints[i]))
            waypoint.extend([0])
            traj.Insert(traj.GetNumWaypoints(), waypoint)
        path.data = [traj]

        joint = path.get_joint(self.robot,0,0)
        self.assertTrue(linalg.norm(joint-self.joints[0])<=1e-5, msg='get_joint failed')

    def test_get_velocity(self):
        """ This test creates a path in openrave format and test the class method get_velocity to return a known joint
        velocity.
        """

        path = blade_coverage.Path(self.rays)
        ind = str()
        for i in range(self.robot.GetDOF()): ind += str(i) + ' '

        cs = ConfigurationSpecification()
        _ = cs.AddGroup('joint_values ' + self.robot.GetName() + ' ' + ind, len(self.robot.GetActiveDOFIndices()),
                        'cubic')
        cs.AddDerivativeGroups(1, False)
        cs.AddDerivativeGroups(2, False)
        _ = cs.AddDeltaTimeGroup()

        traj = RaveCreateTrajectory(self.turbine.env, '')
        traj.Init(cs)
        for i in range(len(self.joints[:1])):
            waypoint = list(self.joints[i])
            waypoint.extend(list(self.joints[i]))
            waypoint.extend(list(self.joints[i]))
            waypoint.extend([0])
            traj.Insert(traj.GetNumWaypoints(), waypoint)
        path.data = [traj]

        velocity = path.get_velocity(self.robot,0,0)
        self.assertTrue(linalg.norm(velocity-self.joints[0])<=1e-5, msg='get_velocity failed')

    def test_get_acc(self):
        """ This test creates a path in openrave format and test the class method get_acc to return a known joint
        acceleration.
        """

        path = blade_coverage.Path(self.rays)
        ind = str()
        for i in range(self.robot.GetDOF()): ind += str(i) + ' '

        cs = ConfigurationSpecification()
        _ = cs.AddGroup('joint_values ' + self.robot.GetName() + ' ' + ind, len(self.robot.GetActiveDOFIndices()),
                        'cubic')
        cs.AddDerivativeGroups(1, False)
        cs.AddDerivativeGroups(2, False)
        _ = cs.AddDeltaTimeGroup()

        traj = RaveCreateTrajectory(self.turbine.env, '')
        traj.Init(cs)
        for i in range(len(self.joints[:1])):
            waypoint = list(self.joints[i])
            waypoint.extend(list(self.joints[i]))
            waypoint.extend(list(self.joints[i]))
            waypoint.extend([0])
            traj.Insert(traj.GetNumWaypoints(), waypoint)
        path.data = [traj]

        acc = path.get_acc(self.robot,0,0)
        self.assertTrue(linalg.norm(acc-self.joints[0])<=1e-5, msg='get_acc failed')

    def test_get_deltatime(self):
        """ This test creates a path in openrave format and test the class method get_deltatime to return a known
        deltatime
        """

        path = blade_coverage.Path(self.rays)
        ind = str()
        for i in range(self.robot.GetDOF()): ind += str(i) + ' '

        cs = ConfigurationSpecification()
        _ = cs.AddGroup('joint_values ' + self.robot.GetName() + ' ' + ind, len(self.robot.GetActiveDOFIndices()),
                        'cubic')
        cs.AddDerivativeGroups(1, False)
        cs.AddDerivativeGroups(2, False)
        _ = cs.AddDeltaTimeGroup()

        traj = RaveCreateTrajectory(self.turbine.env, '')
        traj.Init(cs)
        for i in range(len(self.joints[:1])):
            waypoint = list(self.joints[i])
            waypoint.extend(list(self.joints[i]))
            waypoint.extend(list(self.joints[i]))
            waypoint.extend([0])
            traj.Insert(traj.GetNumWaypoints(), waypoint)
        path.data = [traj]

        deltatime = path.get_deltatime(0,0)
        self.assertTrue(linalg.norm(deltatime)<=1e-5, msg='get_velocity failed')

    def test_get_torques(self):
        """ This test creates a path in openrave format and test the class method get_torques to return a known joint
        torques.
        """

        path = blade_coverage.Path(self.rays)
        ind = str()
        for i in range(self.robot.GetDOF()): ind += str(i) + ' '

        cs = ConfigurationSpecification()
        _ = cs.AddGroup('joint_values ' + self.robot.GetName() + ' ' + ind, len(self.robot.GetActiveDOFIndices()),
                        'cubic')
        cs.AddDerivativeGroups(1, False)
        cs.AddDerivativeGroups(2, False)
        _ = cs.AddDeltaTimeGroup()

        traj = RaveCreateTrajectory(self.turbine.env, '')
        traj.Init(cs)
        for i in range(len(self.joints[:1])):
            waypoint = list(self.joints[i])
            waypoint.extend(list(self.joints[i]))
            waypoint.extend(list(self.joints[i]))
            waypoint.extend([0])
            traj.Insert(traj.GetNumWaypoints(), waypoint)
        path.data = [traj]

        t = [-3.32099384,-77.02582604,123.26938068,-1.00283065,19.89154645,-0.43380835]
        torques = path.get_torques(self.robot,0,0)
        self.assertTrue(linalg.norm(t-torques) <= 1e-5, msg='get_torques failed')

    def test_move_dijkstra(self):
        """ To test move_dijkstra, all new waypoints are verified. i.e., comparison between setUp rays and move_dijkstra
         output.
        """

        path = blade_coverage.Path(self.rays)
        rays = mathtools.equally_spacer(self.rays, self.threshold)
        joint_list, _, _ = path.move_dijkstra(self.turbine, self.rays, self.threshold)
        with self.robot:
            for i in range(len(joint_list)):
                for j in range(len(joint_list[i])):
                    self.robot.SetDOFValues(joint_list[i][j])
                    T = self.robot.GetActiveManipulator().GetTransform()
                    self.assertTrue(linalg.norm(rays[i][j][0:3] - [T[0][3],T[1][3],T[2][3]]) <= 1e-5,
                                    msg='ray verification failed in '+str(i)+','+str(j))
                    n = -array([T[0][0], T[1][0], T[2][0]])
                    n = n/linalg.norm(n)
                    self.assertTrue(dot(rays[i][j][3:6],n)>=cos(self.turbine.config.coating.angle_tolerance),
                                    msg='normal ray verification failed in ' + str(i) + ',' + str(j))

    def test_mls_parallels(self):
        """ To test the mls_parallels an analytical functions q(t) (joints), dq(t) and ddq(t) (derivatives) are computed
        so that |dx| = constant = coating speed.
        dx = J*dq, where J is the Jacobian.
        ddx = J*ddq + dq<sup>T</sup>*H*dq, where H is the Hessian
        dx*ddx = 0 condition for |dx| constant. It can be numerically solved with integral.ode.

        The test will than compare the mls result position, linear velocity and linear accelerations.
        """

        manip = self.robot.GetActiveManipulator()
        times = logspace(0, 1, 100) / 10

        self.robot.SetDOFValues(self.joints[0])
        J = manip.CalculateJacobian()
        v0 = ones(self.robot.GetActiveDOF())*self.turbine.config.coating.coating_speed/linalg.norm(
            dot(J,ones(self.robot.GetActiveDOF())))
        y0, t0 = hstack([self.joints[0], v0]), times[0]

        joints = zeros((len(times), self.robot.GetActiveDOF()))
        joints_velocity = zeros((len(times), self.robot.GetActiveDOF()))
        joints_acc = zeros((len(times), self.robot.GetActiveDOF()))
        joints[0] = y0[:self.robot.GetActiveDOF()]
        joints_velocity[0] = y0[self.robot.GetActiveDOF():]

        def f(t,y):
            with self.robot:
                self.robot.SetDOFValues(y[:self.robot.GetActiveDOF()])
                T = manip.GetTransform()
                J = manip.CalculateJacobian()
                H = self.robot.ComputeHessianTranslation(len(self.robot.GetLinks()) - 1, T[0:3, 3])
            dq = y[self.robot.GetActiveDOF():]
            v = einsum('ij,j,ik->k',J,dq,J)
            v = v/dot(v,v)
            ddq = -v*einsum('ikj,i,j,kl,l',H,dq,dq,J,dq)
            return hstack([dq, ddq])

        r = ode(f).set_integrator('vode', method='bdf')
        r.set_initial_value(y0, t0)

        joints_acc[0] = f(0, y0)[self.robot.GetActiveDOF():]
        for i,t in enumerate(times[1:]):
            r.integrate(t)
            joints[i+1] = r.y[:self.robot.GetActiveDOF()]
            joints_velocity[i+1] = r.y[self.robot.GetActiveDOF():]
            joints_acc[i+1] = f(0,r.y)[self.robot.GetActiveDOF():]

        path = blade_coverage.Path(self.rays)
        MLS_joints, MLS_joint_velocities, MLS_joint_acc = path.mls_parallels(
            self.turbine, [joints], [blade_coverage.compute_dtimes_from_joints(self.turbine,joints)])
        MLS_joints = MLS_joints[0]
        MLS_joint_velocities = MLS_joint_velocities[0]
        MLS_joint_acc = MLS_joint_acc[0]

        for i in range(len(MLS_joints)):
            self.robot.SetDOFValues(joints[i])
            T0 = manip.GetTransform()
            J0 = manip.CalculateJacobian()
            H0 = self.robot.ComputeHessianTranslation(len(self.robot.GetLinks())-1, T0[0:3,3])
            x0 = T0[0:3,3]
            dx0 = dot(J0,joints_velocity[i])
            ddx0 = dot(joints_velocity[i],dot(H0,joints_velocity[i])) + dot(J0,joints_acc[i])

            self.robot.SetDOFValues(MLS_joints[i])
            T1 = manip.GetTransform()
            x1 = T1[0:3,3]
            J1 = manip.CalculateJacobian()
            H1 = self.robot.ComputeHessianTranslation(len(self.robot.GetLinks()) - 1, T1[0:3, 3])
            dx1 = dot(J1, MLS_joint_velocities[i])
            ddx1 = dot(MLS_joint_velocities[i], dot(H1, MLS_joint_velocities[i])) + dot(J1, MLS_joint_acc[i])

            self.assertTrue(linalg.norm(x0-x1) <= 2.5e-2,
                            msg='position verification failed in ' + str(i) + '  ' + str(linalg.norm(x0-x1)) + ' ' + str(x0) + ',' + str(x1))
            self.assertTrue(linalg.norm(dx0-dx1) <= 2.5e-2,
                            msg='joint velocity verification failed in ' + str(i) + '  ' + str(linalg.norm(dx0-dx1)) + ' ' + str(dx0) + ',' + str(dx1) )
            self.assertTrue(linalg.norm(ddx0-ddx1) <= 2.5e-2,
                            msg='joint acc verification failed in ' + str(i) + '  ' + str(linalg.norm(ddx0-ddx1)) + ' ' + str(ddx0) + ',' + str(ddx1))

    def test_mls_joints(self):
        """ This test will validate the error and scale parameters of the MLS, since the mls_joints is a loop to
        decrease the scale in order to reduce the error.
        """

        error = 1e-2
        scale = 3.

        times = cumsum(blade_coverage.compute_dtimes_from_joints(self.turbine, self.joints))/2.
        path = blade_coverage.Path(self.rays)
        new_joints, new_joints_velocities, new_joints_acc = path.mls_joints(
            self.turbine, self.joints, error=error, mls_degree=6, scale=scale, times=times)

        with self.turbine.robot:
            for new_joint, joint in zip(new_joints,self.joints):
                self.turbine.robot.SetDOFValues(joint)
                pos0 = self.turbine.manipulator.GetTransform()[0:3,3]
                self.turbine.robot.SetDOFValues(new_joint)
                pos1 = self.turbine.manipulator.GetTransform()[0:3, 3]
                self.assertTrue(linalg.norm(pos0 - pos1) <= 1e-2,
                                msg='position verification failed: ' + str(pos0) + ',' + str(pos1))


    def test_joint_error_mh12(self):
        """ This is a test specific for the MH12 robot. It will verify if the error between two joints is right"""

        j0 = zeros(self.turbine.robot.GetDOF())
        j1 = j0 + [pi/2, 0, 0, 0, 0, 0]
        p0 = array([1.425, 0., 7.264])
        p1 = array([0., 1.425, 7.264])
        path = blade_coverage.Path(self.rays)
        e, _ = path.joint_error(self.turbine.robot, [j0], [j1])
        error = linalg.norm(p0-p1)
        self.assertAlmostEqual(e, error)

    def test_parallels_transitions(self):

        max_vel = self.turbine.robot.GetDOFMaxVel()
        max_acc = self.turbine.robot.GetDOFMaxAccel()
        dof = self.turbine.robot.GetDOF()
        dof_limits = self.turbine.robot.GetDOFLimits()

        joints_0 = dof_limits[0]/10.
        joints_1 = dof_limits[1]/10.
        vel_0 = max_vel*0.8
        vel_1 = -max_vel*0.8
        acc = zeros(self.turbine.robot.GetDOF())

        joints_parallels = [[joints_0], [joints_1], [joints_0]]
        vel_parallels = [[vel_0], [vel_1], [vel_0]]
        acc_parallels = [[acc], [acc], [acc]]
        times_parallel = [[0], [0], [0]]

        path = blade_coverage.Path(self.rays)
        joints_parallel, joints_vel_parallel, joints_acc_parallel, times_parallel = path.parallels_transitions(
            self.turbine, joints_parallels, vel_parallels, acc_parallels, times_parallel)

        self.assertTrue((abs(vstack(joints_vel_parallel)) <= max_vel*0.8).all())
        #self.assertTrue((abs(vstack(joints_acc_parallel)) <= max_acc * 0.8).all())
        for joints_acc in joints_acc_parallel:
            for acc in joints_acc:
                if (abs(acc) > max_acc * 0.8).any():
                    print abs(acc), max_acc * 0.8

        for i in range(self.turbine.robot.GetDOF()):
            self.assertAlmostEqual(joints_parallel[0][0][i], joints_0[i])
            self.assertAlmostEqual(joints_parallel[2][0][i], joints_1[i])
            self.assertAlmostEqual(joints_parallel[4][0][i], joints_0[i])

        self.assertTrue(len(joints_parallel)==5)

    #def test_get_failed_acc(self):

    #def test_get_failed_vel(self):

    #def test_parallels_transitions(self):

    #def test_create_trajectories(self):

    #def test_organize_rays_in_parallels(self):

    #def test_base_grid_validation(self):

    #def test_compute_dtimes_from_joints(self):

    #def test_compute_dtimes_from_rays(self):


if __name__ == '__main__':
    unittest.main()