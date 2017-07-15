import unittest
from . import TestCase
from .. turbine import Turbine
from .. turbine_config import TurbineConfig
from numpy import zeros, linspace, linalg, array
from .. import blade_coverage
from openravepy import ConfigurationSpecification, RaveCreateTrajectory
import shutil


class TestBladeCoverage(TestCase):
    @classmethod
    def setUpClass(cls):
        super(TestBladeCoverage, cls).setUpClass()
        turbconf = TurbineConfig.load("dummy.cfg", cls.test_dir)
        cls.turbine = Turbine(turbconf)

    def setUp(self):
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

    def test_execute(self):
        path = blade_coverage.Path(self.rays)
        path.execute(self.turbine, 5e-2)
        self.assertTrue(path.success, msg='path was not successfully executed')

    def test_serialize(self):
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
        path = blade_coverage.Path(self.rays)
        try:
            path.deserialize(self.turbine, 'coating_trajectories/test/test_deserialize')
        except:
            self.assertTrue(False, msg='Data can not be deserialized')

    def test_get_joint(self):
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



if __name__ == '__main__':
    unittest.main()