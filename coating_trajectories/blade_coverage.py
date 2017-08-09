from numpy import array, linalg, dot, zeros, inf, vstack, mean, std, cumsum, abs, sign, ones, cross, linspace
from numpy.polynomial import legendre
from math import acos
import planning
from openravepy import ConfigurationSpecification, interfaces, planningutils, RaveCreateTrajectory, interfaces
import mathtools
import errno
import rail_place
from xml.etree import ElementTree as ET
from os import listdir, makedirs
from os.path import realpath, join, isfile
import time as Time
import matplotlib.pyplot as plt


## @file
# @brief This contains functions and a class (path) to compute the joint solutions given trajectories (operational to joint space)
# @author Renan S. Freitas
# @bug No known bugs

class Path:
    """ Class to compute, store, and serialize joint path.
        Robot's full trajectories are stored in OpenRave format Trajectory Class, as below:
        [joints_values, joints_vel_values, joints_acc_values, deltatimes]

        Args:
            rays: (float[m][n<SUB>i</SUB>][6]) cartesian points and normals

        Examples
            >>> res = path(rays)
    """

    def __init__(self, rays=None):
        self.rays = rays
        self.data = []
        self.success = False

    def execute(self, turbine, threshold=5e-2, dtimes=None):
        """ Method to compute joint_values, joint_velocities, joint_accelerations and deltatimes.
            It uses the Dijkstra planning algorithm (see move_dijkstra function).

        Args:
            turbine: (@ref Turbine) is the turbine object.
            blade: (@ref Turbine) is the blade object (you might use DB.load_blade())
            threshold: (float) is the interpolation threshold, as rays are usually well spaced.

        Examples:
            >>> path.execute(turbine, DB.load_blade())
        """

        if self.rays is None: return

        joint_path, rays_list, dtimes = self.move_dijkstra(turbine, self.rays, threshold)
        if len(joint_path) == 0:
            return

        joint_path, dtimes = self.parallel_transition(turbine, joint_path, dtimes)

        new_joint_path, new_joint_velocity_path, new_joint_acc_path, new_times_path = self.mls_parallels(
            turbine, joint_path, dtimes)

        acc = abs(vstack(new_joint_acc_path))
        vel = abs(vstack(new_joint_velocity_path))

        if (vel > turbine.robot.GetDOFMaxVel()).any():
            print 'vel max fail'
            return

        if (acc > turbine.robot.GetDOFMaxAccel()).any():
            print 'acc max fail'
            return

        self.success = True

        self.data = []
        for i in range(len(new_joint_path)):
            traj = self.create_trajectory(turbine,
                                          new_joint_path[i],
                                          new_joint_velocity_path[i],
                                          new_joint_acc_path[i],
                                          new_times_path[i])
            self.data.append(traj)

        return

    def serialize(self, directory=''):
        """ Method serializes data in OpenRave format output. It is an xml file.

        Args:
            directory: (str) is the relative path to the folder where to save.

        Examples:
            >>> path.serialize('my_folder')
        """

        try:
            makedirs(directory)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

        for i, traj in enumerate(self.data):
            tree = ET.XML(traj.serialize())
            with open(join(directory, str(i)), 'w') as f:
                f.write(ET.tostring(tree))
        return

    def deserialize(self, turbine, directory):
        """ Method deserializes data in OpenRave format.

        Args:
            turbine: (@ref Turbine) is the turbine object.
            directory: (str) is the relative path to the folder where to load.

        Examples:
            >>> path.deserialize(turbine,'my_folder')
        """

        ind = str()
        for i in range(turbine.robot.GetDOF()): ind += str(i) + ' '

        cs = ConfigurationSpecification()
        _ = cs.AddGroup('joint_values ' + turbine.robot.GetName() + ' ' + ind, len(turbine.robot.GetActiveDOFIndices()),
                        'cubic')
        cs.AddDerivativeGroups(1, False)
        cs.AddDerivativeGroups(2, False)
        _ = cs.AddDeltaTimeGroup()

        directory = realpath(directory)
        TRAJ = []
        onlyfiles = [f for f in listdir(directory) if isfile(join(directory, f))]
        onlyfiles.sort(key=int)
        for afile in onlyfiles:
            traj = RaveCreateTrajectory(turbine.env, '')
            traj.Init(cs)
            xml = ET.parse(open(join(directory, afile)))
            traj.deserialize(ET.tostring(xml.getroot()))
            TRAJ.append(traj)
        self.data = TRAJ
        return

    def simulate(self, robot, parallel_number):
        """ Method simulates and call visualization, in real time: the robot performing coat.

        Args:
            robot: (Turbine.robot) is the robot object.
            parallel_number: (int) is the parallel number to be simulated.

        Examples:
            >>> path.simulate(turbine.robot, 0)
            >>> for i in range(len(path.data)): path.simulate(turbine.robot, i)
        """

        _ = robot.GetController().SetPath(self.data[parallel_number])
        robot.WaitForController(0)
        return

    def get_joint(self, robot, parallel_number, point_number):
        """ Method gets a specific joint_value from path.

        Args:
            robot: (Turbine.robot) is the robot object.
            parallel_number: (int) is the parallel number to get the joint.
            point_number: (int) is the specific point in the parallel to get the joint _value.

        Returns:
            List float[m][n<SUB>i</SUB>][nDOF] of joint values

        Examples:
            >>> path.get_joint(turbine.robot, 0, 0)
            >>> N = path.data[0].GetNumWaypoints()
            >>> joints = []
            >>> for i in range(N): joints.append(path.get_joint(turbine.robot,0,i))
        """

        traj = self.data[parallel_number]
        spec = traj.GetConfigurationSpecification()
        return spec.ExtractJointValues(traj.GetWaypoint(point_number), robot, range(robot.GetDOF()))

    def get_velocity(self, robot, parallel_number, point_number):
        """ Method gets a specific joint_velocity from path.

        Args:
            robot: (Turbine.robot) is the robot object.
            parallel_number: (int) is the parallel number to get the joint.
            point_number: (int) is the specific point in the parallel to get the joint _value.

        Returns:
            List float[m][n<SUB>i</SUB>][nDOF] of joint values

        Examples:
            >>> path.get_velocity(turbine.robot, 0, 0)
            >>> N = path.data[0].GetNumWaypoints()
            >>> velocities = []
            >>> for i in range(N): velocities.append(path.get_velocity(turbine.robot,0,i))
        """

        traj = self.data[parallel_number]
        spec = traj.GetConfigurationSpecification()
        return spec.ExtractJointValues(traj.GetWaypoint(point_number), robot, range(robot.GetDOF()), 1)

    def get_acc(self, robot, parallel_number, point_number):
        """ Method gets a specific joint_acceleration from path.

        Args:
            robot: (Turbine.robot) is the robot object.
            parallel_number: (int) is the parallel number to get the joint.
            point_number: (int) is the specific point in the parallel to get the joint _value.

        Returns:
            List float[m][n<SUB>i</SUB>][nDOF] of joint values

        Examples:
            >>> path.get_acc(turbine.robot, 0, 0)
            >>> N = path.data[0].GetNumWaypoints()
            >>> accelerations = []
            >>> for i in range(N): accelerations.append(path.get_acc(turbine.robot,0,i))
        """

        traj = self.data[parallel_number]
        spec = traj.GetConfigurationSpecification()
        return spec.ExtractJointValues(traj.GetWaypoint(point_number), robot, range(robot.GetDOF()), 2)

    def get_deltatime(self, parallel_number, point_number):
        """ Method gets a specific joint_acceleration from path.

        Args:
            parallel_number: (int) is the parallel number to get the joint.
            point_number: (int) is the specific point in the parallel to get the joint _value.

        Returns:
            List float[m][n<SUB>i</SUB>][nDOF] of joint values

        Examples:
            >>> path.get_deltatime(0, 0)
            >>> N = path.data[0].GetNumWaypoints()
            >>> times = []
            >>> for i in range(N): times.append(path.get_deltatime(0,i))
        """

        traj = self.data[parallel_number]
        spec = traj.GetConfigurationSpecification()
        return spec.ExtractDeltaTime(traj.GetWaypoint(point_number))

    def get_torques(self, robot, parallel_number, point_number):
        """ Method gets a specific joint_acceleration from path.

            Args:
                robot: (Robot) is the robot object.
                parallel_number: (int) is the parallel number to get the joint.
                point_number: (int) is the specific point in the parallel to get the joint _value.

            Returns:
                List float[m][n<SUB>i</SUB>][nDOF] of joint values

            Examples:
                >>> path.get_torques(turbine.robot, 0, 0)
                >>> N = path.data[0].GetNumWaypoints()
                >>> torques = []
                >>> for i in range(N): torques.append(path.get_torques(turbine.robot,0,i)) # Get all torques in parallel 0
        """

        with robot:
            robot.SetDOFValues(self.get_joint(robot, parallel_number, point_number))
            robot.SetDOFVelocities(self.get_velocity(robot, parallel_number, point_number))
            return robot.ComputeInverseDynamics(self.get_acc(robot, parallel_number, point_number))

    def plot_velocities(self, turbine, parallel_number):
        velocities = []
        dtimes = []
        max_vel = turbine.robot.GetDOFMaxVel()

        N = self.data[parallel_number].GetNumWaypoints()

        for i in range(N):
            velocities.append(self.get_velocity(turbine.robot, parallel_number, i))
            dtimes.append(self.get_deltatime(parallel_number, i))

        velocities = array(velocities)
        f, ax = plt.subplots(turbine.robot.GetDOF(), sharex=True)
        dtimes = cumsum(dtimes)

        for i in range(turbine.robot.GetDOF()):
            ax[i].plot(dtimes, velocities[:, i], color='b')
            ax[i].plot(dtimes, ones(N) * max_vel[i], color='r')
            ax[i].plot(dtimes, -ones(N) * max_vel[i], color='r')
            ax[i].set_title('Joint Velocity ' + str(i))

        f.subplots_adjust(hspace=0.3)
        plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
        plt.show()
        return

    def plot_acc(self, turbine, parallel_number):
        acc = []
        dtimes = []
        max_acc = turbine.robot.GetDOFMaxAccel()

        N = self.data[parallel_number].GetNumWaypoints()

        for i in range(N):
            acc.append(self.get_acc(turbine.robot, parallel_number, i))
            dtimes.append(self.get_deltatime(parallel_number, i))

        acc = array(acc)
        f, ax = plt.subplots(turbine.robot.GetDOF(), sharex=True)
        dtimes = cumsum(dtimes)

        for i in range(turbine.robot.GetDOF()):
            ax[i].plot(dtimes, acc[:, i], color='b')
            ax[i].plot(dtimes, ones(N) * max_acc[i], color='r')
            ax[i].plot(dtimes, -ones(N) * max_acc[i], color='r')
            ax[i].set_title('Joint Acc ' + str(i))

        f.subplots_adjust(hspace=0.3)
        plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
        plt.show()
        return

    @staticmethod
    def move_dijkstra(turbine, organized_rays_list, interpolation):
        """ Given parallels (cartesian points x-y-z-nx-ny-nz), this method returns robot's joints solution, using
        inverse kinematics (openrave ikfast) and the Dijktra algorithm for path optimization:
        @see https://en.wikipedia.org/wiki/Dijkstra%27s_algorithm

        Args:
            turbine: (@ref Turbine) turbine object
            organized_rays_list: (float[m][n<SUB>i</SUB>][6]) zigzagging parallels
            interpolation: (float) distance between points in parallels

        Returns:
             Two lists: joint values and equally spaced zigzagging parallels

        Examples:
            >>> joint_path, rays_list = path.move_dijkstra(turbine, rays_list, 3e-3)
        """

        robot = turbine.robot
        joint_path_list = []
        vel_limits = array(robot.GetDOFVelocityLimits())[:-1]
        acc_limits = array(robot.GetDOFMaxAccel())[:-1]
        organized_rays_list = mathtools.equally_spacer(organized_rays_list, interpolation)[:3]
        dtimes_list = []

        for i, organized_rays in enumerate(organized_rays_list):

            if i != len(organized_rays_list) - 1:
                organized_rays += [organized_rays_list[i + 1][0]]

            try:
                joints = planning.joint_planning(turbine, organized_rays, True)
            except IndexError:
                continue

            if i != 0:
                joints.insert(0, [joint_path_list[-1][-1]])

            dtimes = compute_dtimes_from_joints(turbine, [j[0] for j in joints])
            joints = array([array(j)[:, :-1] for j in joints])
            joint_path, path, min_cost, adj, cost = planning.make_dijkstra_vel(joints, dtimes, vel_limits,
                                                                           acc_limits, True)
            print  min_cost

            if i != len(organized_rays_list) - 1:
                joint_path = joint_path[:-1]
                dtimes = dtimes[:-1]

            if i != 0:
                joint_path = joint_path[1:]
                dtimes = dtimes[1:]

            dtimes_list.append(dtimes)
            joint_path_complete = zeros((len(joint_path), robot.GetDOF()))
            joint_path_complete[:len(joint_path), :len(joint_path[0])] = joint_path
            joint_path_list.append(joint_path_complete)

        return joint_path_list, organized_rays_list, dtimes_list

    @staticmethod
    def refine_dijkstra(turbine, joint_path_list, rays_list, interpolation):
        """ Method to refine the dijkstra algorithm. Since there is orientation tolerance for the end-effector,
        there might be infinite inverse kinematics solutions for all points, thus refine dijkstra will try to find
        closer (optimal) solutions, using the previous dijkstra general solution (rough discretization) .

        Args:
            turbine: (@ref Turbine) turbine object
            joint_path_list: (float[m][n<SUB>i</SUB>][nDOF]) joint values list for each parallel
            rays_list: (float[m][n<SUB>i</SUB>][6]) zigzagging equally spaced parallels
            interpolation: (float) distance between points

        Returns:
            New joint values list for each parallel

        Examples:
            >>> joint_path = path.refine_dijkstra(turbine, joint_path_list, rays_list, 3e-3)
        """
        robot = turbine.robot
        time = interpolation / turbine.config.coating.coating_speed
        limits = robot.GetDOFVelocityLimits() * time
        deep = True
        new_joint_path = []
        with robot:
            for i, rays in enumerate(rays_list):
                joints_path = joint_path_list[i]
                new_joints = []
                for j, joint in enumerate(joints_path):
                    if j == 0 or j == len(joints_path) - 1:
                        new_joints.append([joint])
                        continue

                    robot.SetDOFValues(joints_path[j - 1])
                    Rx = robot.GetActiveManipulator().GetTransform()[0:3, 0]
                    Rx = Rx / linalg.norm(Rx)
                    d = -dot(rays[j - 1][3:6], Rx)
                    angle0 = mathtools.clean_acos(d)

                    robot.SetDOFValues(joint)
                    Rx = robot.GetActiveManipulator().GetTransform()[0:3, 0]
                    Rx = Rx / linalg.norm(Rx)
                    d = -dot(rays[j][3:6], Rx)
                    angle1 = mathtools.clean_acos(d)

                    angle_tolerance_init = max([min([angle0, angle1]) - 0.01, 0])
                    angle_tolerance_end = min([max([angle0, angle1]) + 0.01, 3.14])
                    iksol = planning.ik_angle_tolerance(turbine, rays[j],
                                                        angle_tolerance_init=angle_tolerance_init,
                                                        angle_tolerance_end=angle_tolerance_end,
                                                        number_of_phi=36, number_of_theta=8, deep=deep)
                    iksol += [joint]
                    new_joints.append(iksol)

                joint_path, path, min_cost, adj, cost = planning.make_dijkstra_mh12(new_joints, limits, True)
                if min_cost != inf:
                    new_joint_path.append(joint_path)
                else:
                    return joint_path_list
        return new_joint_path

    def mls_parallels(self, turbine, joints_parallels, dtimes_parallels):
        """ This method calls mls_joints to smooth the trajectories (all parallels)

        Args:
            turbine: (@ref Turbine) turbine object
            joint_path: (float[m][n<SUB>i</SUB>][nDOF]) list of joints for all parallels

        Returns:
            Four lists: smooth joint values, joint velocities,  joint accelerations, and deltatimes

        Examples:
            >>> joint_path, joint_velocity_path, joint_acc_path, dtimes = path.mls_parallels(turbine, joint_path, dtimes_path)
        """

        new_joints_parallels = []
        new_joints_vel_parallels = []
        new_joints_acc_parallels = []
        new_dtimes_parallels = []

        for joints, dtimes in zip(joints_parallels, dtimes_parallels):
            scale = 3.
            dtimes = list(dtimes)
            self.mls_joints(turbine, joints, scale, dtimes)

            joints = array( (2*joints[0]-joints[1]) + list(joints) + 2*(joints[-1]-joints[-2]) )
            dtimes = dtimes[:2] + dtimes[1:] + dtimes[-1]
            new_joints, new_joints_velocities, new_joints_acc, new_dtimes = self.mls_joints(
                turbine, joints, scale, dtimes)

            new_joints_parallels += [new_joints]
            new_dtimes_parallels += [new_dtimes]
            new_joints_vel_parallels += [new_joints_velocities]
            new_joints_acc_parallels += [new_joints_acc]

        return new_joints_parallels, new_joints_vel_parallels, new_joints_acc_parallels, new_dtimes_parallels

    def mls_joints(self, turbine, joints, error=2.5e-3, mls_degree = 6, scale=3., times=None):
        """ The Dijkstra optimization method will calculate the best path for the given discretization. Discretization
        will generate some "jumps", and infeasible velocities/accelerations. A moving least square method was developed
        for path smoothness. The MLS generates a trajectory that will not pass by all waypoints, but it will choose the
        best scale given the maximum required error in cartesian space (meters).
        It returns joints, joint velocities, accelerations and deltatimes.
        This method chooses the best scale given the maximum required error in cartesian space.

        Args:
            turbine: (@ref Turbine) turbine object
            joints: (float[n<SUB>i</SUB>][nDOF]) list of joints]
            error: (float) maximum required error in cartesian space (meters)
            scale: (float) MLS scale
            times: (float[n<SUB>i</SUB>]) deltatimes between joints

        Returns:
            lists: smooth joint values, joint velocities,  joint accelerations, and deltatimes

        Examples:
            >>> joints, joint_velocity, joint_acc, dtimes = path.mls_joints(turbine, joints)
        """

        if times is None:
            times = cumsum(compute_dtimes_from_joints(turbine, joints))

        while scale > 0:
            new_joints = []
            new_joints_velocities = []
            new_joints_acc = []
            for i, joint in enumerate(array(joints).T):
                j, v, a = mathtools.legMLS(joint, times, times, mls_degree, scale)
                new_joints += [j]
                new_joints_velocities += [v]
                new_joints_acc += [a]
            new_joints = array(new_joints).T
            new_dtimes = compute_dtimes_from_joints(turbine, new_joints)
            new_joints_velocities = array(new_joints_velocities).T
            new_joints_acc = array(new_joints_acc).T
            e, _ = self.joint_error(turbine.robot, joints, new_joints)
            if e <= error:
                break
            scale -= .1
        return array(new_joints), array(new_joints_velocities), array(new_joints_acc), array(new_dtimes)

    @staticmethod
    def joint_error(robot, joints_a, joints_b):
        """ This method calculates the mean(error) + 0.5*std(error) between two lists of joints (a and b) in cartesian
        space. It iterates in lists joints_b and joints_a, setting the robot DOF values and comparing the end-effector
        position. Orientations are not compared.

        Args:
            robot: (Robot) object
            joints_a: (float[n<SUB>i</SUB>][nDOF]) list of joints_a
            joints_b: (float[n<SUB>i</SUB>][nDOF]) list of joints_b to be compared with.

        Returns:
            mean(error) + 0.5*std (float).

        Examples:
            >>> e = path.joint_error(robot, joints_a, joints_b)
        """

        manip = robot.GetActiveManipulator()
        error = []
        new_points = []
        with robot:
            for j in range(len(joints_b)):
                robot.SetDOFValues(joints_b[j])
                P0 = manip.GetTransform()[0:3, 3]
                robot.SetDOFValues(joints_a[j])
                P1 = manip.GetTransform()[0:3, 3]
                error.append(linalg.norm(P0 - P1))
                new_points += [P0]
        return mean(error) + .5 * std(error), new_points

    def replanning(self, turbine, new_joint_path, new_joint_velocity_path, new_joint_acc_path, new_dtimes_path):
        """ Not implemented.
        """

        for k in range(len(new_joint_path)):
            q, dq, ddq, dt = array(new_joint_path[k]), array(new_joint_velocity_path[k]), array(
                new_joint_acc_path[k]), array(new_dtimes_path[k])

            for j in range(len(ddq[0])):
                timeshift_list = []
                # start_h = -1
                # start = -1
                max_accel = turbine.robot.GetDOFMaxAccel()[j]
                for i, ddqij in enumerate(ddq[:, j]):
                    #     if ddqij >  max_accel*0.95:
                    #         if start == -1:
                    #             start_h = i
                    #     if ddqij >  max_accel:
                    #         if start == -1:
                    #             start = 1
                    #     if ddqij <  max_accel*0.95:
                    #         if start == 1:
                    #             timeshift_list.append([start_h, i])
                    #             start = -1
                    #         start_h = -1

                    if abs(ddqij) > max_accel * 0.95:
                        timeshift_list.append(i)

                if len(timeshift_list) > 0:
                    q[:, j], dt = self.timeshift(timeshift_list, max_accel * 0.95, q[:, j], dq[:, j], ddq[:, j], dt)
                    try:
                        q, dq, ddq, dt = self.mls_joints(turbine, q, 3., cumsum(dt))
                    except TypeError:
                        print 'dt type - ', type(dt)

            new_joint_path[k], new_joint_velocity_path[k], new_joint_acc_path[k], new_dtimes_path[k] = q, dq, ddq, dt
        return new_joint_path, new_joint_velocity_path, new_joint_acc_path, new_dtimes_path

    @staticmethod
    def timeshift(timeshift_list, max_accel, q, dq, ddq, dt):
        """ Not implemented.
        """

        # for times in timeshift_list:
        #     dq0 = dq[times[0]]
        #     delta_q = q[times[1]]-q[times[0]]
        #     dtime = sum(dt[times[0] + 1:times[1] + 1])
        #     limit = sign(ddq[times[0]]) * min(max_accel,abs(dq[times[1]]-dq[times[0]])/dtime)
        #     tf = (-dq0 + (dq0 ** 2 + 2 * delta_q * limit) ** 0.5) / limit
        #
        #     print 'times[0],times[1] - ', times[0], times[1]
        #     print 'ddq[times[0]] - ', ddq[times[0]]
        #     print 'ddq[times[1]] - ', ddq[times[1]]
        #     print 'delta_dq/dtime - ', (dq[times[1]]-dq[times[0]])/dtime
        #     print 'dq0 - ', dq0
        #     print 'delta_q - ', delta_q
        #     print 'limit - ', limit
        #     print 'dq0**2+2*delta_q*limit - ', dq0 ** 2 + 2 * delta_q * limit
        #     print '(dq0**2+2*delta_q*limit)**0.5 -  ', (dq0 ** 2 + 2 * delta_q * limit) ** 0.5
        #     print '(-dq0+(dq0**2+2*delta_q*limit)**0.5) - ', (-dq0 + (dq0 ** 2 + 2 * delta_q * limit) ** 0.5)
        #
        #     dt[times[0]:times[1]] = tf / (times[1] - times[0])
        #     q[times[0]:times[1]] = linspace(q[times[0]],q[times[1]],times[1]-times[0],endpoint=False)

        for t in timeshift_list:
            if t == len(q) - 1:
                continue

            dq0 = dq[t]
            delta_q = q[t + 1] - q[t]
            limit = sign(ddq[t]) * max_accel * 0.95
            tf = (-dq0 + (dq0 ** 2 + 2 * delta_q * limit) ** 0.5) / limit

            print 't - ', t
            print 'ddq[t] - ', ddq[t]
            print 'limit - ', limit
            print 'dq0 - ', dq0
            print 'dqf - ', dq[t + 1]
            print 'delta_q - ', delta_q
            print 'tf - ', tf
            print 'dt[t+1] -', dt[t + 1]

            dt[t + 1] = tf

        return q, dt

    @staticmethod
    def parallels_transitions(turbine, joints_parallel, times_parallel, step=0.1, number_of_points=6., max_acc=.8) :
        """ Given all parallels solutions, i.e., joints solutions for all waypoints split in parallels (lists),
        this method computes cubic polynomials to interpolate the two points: end of a parallel - begin of the next
        parallel, in joint space, considering velocities. It will return the parallels and the computed transition
        paths between them.

        Args:
            turbine: (@ref Turbine) turbine object
            joints_parallel: (float[m][n<SUB>i</SUB>][nDOF]) list of joints for all parallels
            times_parallel: (float[m][n<SUB>i</SUB>]) list of deltatimes for all parallels
            step: distance (meters) between points
            number_of_points: minimal number of points
            max_acc: maximum permitted acceleration (percentage)

        Returns:
            Modifies the current joints_parallel and times_parallel with the computed transitions.
        """

        robot = turbine.robot

        for i in range(len(joints_parallel) - 1):
            joint_0 = joints_parallel[2*i][-1]
            vel_0 = (joints_parallel[2*i][-1]-joints_parallel[i][-2])/times_parallel[2*i][-1]
            joint_1 = joints_parallel[2*i + 1][0]
            joint_2 = joints_parallel[2*i + 1][1]
            dt12 = times_parallel[2*i + 1][1]

            dq = (joint_2 - joint_1)/dt12

            tf = max((dq - vel_0) / (robot.GetDOFMaxAccel() * max_acc))

            c = mathtools.legcubic_path(joint_0, joint_1, vel_0, dq, [0, tf])

            joints = []
            times = []
            dt = min([step, tf / number_of_points])
            for t in linspace(0, tf, max([tf / step, number_of_points])):
                joints.append(legendre.legval(t, c))
                times.append(dt)

            joints = joints[1:]
            joints_parallel[2 * i + 1] =  joints_parallel[2*i+1][1:]
            times = times[1:]
            times_parallel[2 * i + 1] = times_parallel[2*i+1][1:]
            joints_parallel.insert(2*i + 1, joints)
            times_parallel.insert(2*i + 1, times)

        return joints_parallel, times_parallel

    @staticmethod
    def create_trajectory(turbine, joints, joints_vel, joints_acc, times):
        """ This method creates the trajectory specification in OpenRave format and insert the waypoints:

        joints - vels - accs - times
        DOF - DOF - DOF - DOF - 1

        Args:
            turbine: (@ref Turbine) turbine object
            joints:  (float[n<SUB>i</SUB>][nDOF]) joints to create the trajectory
            joints_vel: (float[n<SUB>i</SUB>][nDOF]) joints_vel to create the trajectory
            joints_acc: (float[n<SUB>i</SUB>][nDOF]) joints_acc to create the trajectory
            times: (float[n<SUB>i</SUB>]) deltatimes

        Returns:
            trajectory: OpenRave object.
        """

        robot = turbine.robot

        ind = str()
        for i in range(robot.GetDOF()): ind += str(i) + ' '

        cs = ConfigurationSpecification()
        _ = cs.AddGroup('joint_values ' + robot.GetName() + ' ' + ind,
                        len(robot.GetActiveDOFIndices()),
                        'cubic')
        cs.AddDerivativeGroups(1, False)
        cs.AddDerivativeGroups(2, False)
        _ = cs.AddDeltaTimeGroup()

        traj = RaveCreateTrajectory(turbine.env, '')
        traj.Init(cs)

        for i in range(len(joints)):
            waypoint = list(joints[i])
            waypoint.extend(list(joints_vel[i]))
            waypoint.extend(list(joints_acc[i]))
            waypoint.extend([times[i]])
            traj.Insert(traj.GetNumWaypoints(), waypoint)
        return traj

def organize_rays_in_parallels(DB, grid):
    """ Function makes a zigzagging path from parallels.

    Args:
        DB: (@ref DB) database object.
        grid: (int) grid to be coated.

    Returns:
        float[m][n<SUB>i</SUB>][6] organized rays to be coated (zigzagging).

    Examples:
        >>> organized_rays = organize_rays_in_parallels(DB, 0)
    """

    rays = DB.compute_rays_from_grid(grid)

    organized_rays = []
    for i in range(0, len(rays)):
        not_empty_rays = mathtools.notempty(rays[i])
        if len(not_empty_rays) == 0:
            continue
        if i % 2 == 0:
            organized_rays.append(not_empty_rays)
        else:
            organized_rays.append(list(reversed(not_empty_rays)))
    return organized_rays


def base_grid_validation(turbine, psa, DB, grid, threshold=5e-2):
    """ Given the real blade angle, this method tries to validate a base:
    1) rotate the blades (update the environment);
    2) organize trajectories (removing empty points, adding borders,
    and making zigzagging lists);
    3) Dijkstra algorithm.
    4) Moving Least Square smoothness;

    Args:
        turbine: (@ref Turbine) is the turbine object.
        psa: (tuple[3]) primary, secondary, alpha (base position)
        DB: (@ref DB) database object.
        grid: (int) grid to be coated.
        threshold: (float) is the interpolation threshold, as rays are usually well spaced.

    Returns:
        path object

    Examples:
        >>> path = base_grid_validation(turbine, psa, DB, 0)
    """

    T = DB.T
    T = dot(T, linalg.inv(turbine.blades[3].GetTransform()))  # HARDCODED

    for blade in turbine.blades:
        T0 = blade.GetTransform()
        blade.SetTransform(dot(T, T0))

    rp = rail_place.RailPlace(psa)
    turbine.place_rail(rp)
    turbine.place_robot(rp)
    organized_rays_list = organize_rays_in_parallels(DB, grid)
    trajectory = Path(organized_rays_list)
    trajectory.execute(turbine, threshold)
    return trajectory


def compute_dtimes_from_joints(turbine, joints):
    """ Given the joints solution, this method set robot DOF values and access the config file (velocity requirements),
    computing the required delta times between joints.

    Args:
        turbine: (@ref Turbine) is the turbine object.
        joints: (float[n][nDOF]) robot joints for a parallel

    Returns:
        a list of delta times

    Examples:
        >>> dtimes = compute_dtimes_from_joints(turbine, joints)
    """

    rays = []

    with turbine.robot:
        for joint in joints:
            turbine.robot.SetDOFValues(joint)
            rays += [turbine.robot.GetActiveManipulator().GetTransform()[:3, 3]]

    return compute_dtimes_from_rays(turbine, rays)


def compute_dtimes_from_rays(turbine, rays):
    """ Given the rays, this method access the config file (velocity requirements) and compute the required
        delta times between joints.

    Args:
        turbine: (@ref Turbine) is the turbine object.
        rays: (float[n][6]) x,y,z,nx,ny,nz cartesian coordinates and normal vector

    Returns:
        a list of delta times

    Examples:
        >>> dtimes = compute_dtimes_from_joints(turbine, joints)
    """
    rays = array(rays)
    v = turbine.config.coating.coating_speed

    return array([0.] + list(linalg.norm(rays[1:, 0:3] - rays[:-1, 0:3], axis=1) / v))