from numpy import array, linalg, dot, zeros, vstack, mean, std, cumsum, abs, ones, linspace, clip, vander, arange
from numpy.polynomial import legendre
from copy import copy

import dijkstra_vel
import dijkstra_acc
import robot_utils
from openravepy import ConfigurationSpecification, interfaces, planningutils, RaveCreateTrajectory, interfaces
import mathtools
import errno
import rail_place
from xml.etree import ElementTree as ET
from os import listdir, makedirs
from os.path import realpath, join, isfile
import matplotlib.pyplot as plt
from time import sleep


## @file
# @brief This contains functions and a class (path) to compute the joint solutions given trajectories (operational to joint space)
# @author Renan S. Freitas & Eduardo Elael
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
        self.joint_dijkstra = []
        self.dtimes = []
        self.iksols = []

    def execute(self, turbine, threshold=5e-2):
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

        self.compute_iksol(turbine, self.rays, threshold)

        self.joint_dijkstra, self.dtimes = self.move_dijkstra(turbine, self.iksols)

        joint_path, joint_velocity, joint_acc = self.mls_parallels(
            turbine, self.joint_dijkstra, self.dtimes)

        self.create_trajectory_data(turbine, joint_path, joint_velocity, joint_acc)

        #new_joint_path, new_joint_velocity_path, new_joint_acc_path, dtimes = self.parallels_transitions(
        #    turbine, new_joint_path, new_joint_velocity_path, new_joint_acc_path, dtimes)


        acc_middle = []
        for acc_path in joint_acc:
            acc_middle.append(acc_path[1:-1])

        vel_middle = []
        for vel_path in joint_velocity:
            vel_middle.append(vel_path[1:-1])

        acc = abs(vstack(acc_middle))
        vel = abs(vstack(vel_middle))

        if (vel > turbine.robot.GetDOFMaxVel()).any():
            self.vel_max_fail = True

        if (acc > turbine.robot.GetDOFMaxAccel()).any():
            self.acc_max_fail = True
            # A further inspection must be made. There are some strategies: break the parallel, retiming

        self.grid_break, self.success = self.break_grid(turbine)
        return


    def create_trajectory_data(self, turbine, pos, velocity, acc):
        self.data = []
        for i in range(len(pos)):
            traj = self.create_trajectory(turbine,
                                          pos[i],
                                          velocity[i],
                                          acc[i],
                                          self.dtimes[i])
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

    def simulate(self, robot, vis, n=0):
        """ Method simulates and call visualization, in real time: the robot performing coat.
        Args:
            robot: (Turbine.robot) is the robot object.
            parallel_number: (int) is the parallel number to be simulated.
            n: (int) interpolate n points between data
        Examples:
            >>> path.simulate(turbine.robot, 0)
            >>> for i in range(len(path.data)): path.simulate(turbine.robot, i)
        """

        manip = robot.GetActiveManipulator()
        if n == 0:
            for parallel_number in range(len(self.data)):
                for point_number in range(self.data[parallel_number].GetNumWaypoints()):
                    robot.SetDOFValues(self.get_joint(robot, parallel_number, point_number))
                    p = manip.GetTransform()[0:3,3]
                    _ = vis.plot(p, 'point', (1,0,0))
                    sleep(self.get_deltatime(parallel_number, point_number))

        else:
            init = self.get_joint(robot, 0, 0)
            robot.SetDOFValues(init)
            for parallel_number in range(len(self.data)):
                for point_number in range(self.data[parallel_number].GetNumWaypoints()):
                    end = self.get_joint(robot, parallel_number, point_number)
                    dtime = self.get_deltatime(parallel_number, point_number)
                    if point_number > 0:
                        joints = mathtools.linspace_array(init, end, n)
                        for joint in joints:
                            robot.SetDOFValues(joint)
                            p = manip.GetTransform()[0:3, 3]
                            _ = vis.plot(p, 'point', (0, 0, 1),10)
                            sleep(dtime/n)
                    else:
                        robot.SetDOFValues(end)
                        p = manip.GetTransform()[0:3, 3]
                        _ = vis.plot(p, 'point', (0, 0, 1),10)
                        sleep(dtime)
                    init = end


        #_ = robot.GetController().SetPath(self.data[parallel_number])
        #robot.WaitForController(0)
        return

    def plot_coated_points(self, robot, vis, step = 10):
        points = []
        with robot:
            for parallel_index in range(len(self.data)):
                for point_index in range(self.data[parallel_index].GetNumWaypoints()):
                    joint1 = self.get_joint(robot, parallel_index, point_index)
                    robot.SetDOFValues(joint1)
                    points.append(robot.GetActiveManipulator().GetTransform()[0:3,3])

        for i in range(len(points)-1):
            x = linspace(points[i][0],points[i+1][0], step)
            y = linspace(points[i][1],points[i+1][1], step)
            z = linspace(points[i][2],points[i+1][2], step)
            points.append(zip(x,y,z))

        _ = vis.plot_lists(points,'interpolation',(1,0,0))
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

    def get_linear_acc(self, robot, parallel_number, point_number, link):
        """Method returns end-effector acc.

        Args:
            robot: (Robot) is the robot object.
            parallel_number: (int) is the parallel number to get the joint.
            point_number: (int) is the specific point in the parallel to get the joint _value.
            link: link to get linear acc. (-3 for last link in mh12 with flame and gun)

        Examples:
            >>> path.get_linear_acc(turbine.robot, 0, 0)
            >>> N = path.data[0].GetNumWaypoints()
            >>> lin_acc = []
            >>> for i in range(N): lin_acc.append(path.get_linear_acc(turbine.robot,0,i)) # Get all linear acc in parallel 0
        """

        with robot:
            robot.SetDOFValues(self.get_joint(robot, parallel_number, point_number))
            robot.SetDOFVelocities(self.get_velocity(robot, parallel_number, point_number))
            return robot.GetLinkAccelerations(self.get_acc(robot,parallel_number,point_number))[link]



    def get_failed_acc(self, robot):
        """ Method returns a list [parallel_number, point_number] (Path.data structure) with all accelerations above
        the limit.
        Args:
            robot: (Robot) is the robot object.
        Returns:
            List int[m][2] [parallel_number, point_number].
        Examples:
            >>> acc_fail = path.get_failed_acc(robot)
        """

        acc_not_valid = []
        for parallel_number in range(len(self.data)):
            for point_number in range(1,self.data[parallel_number].GetNumWaypoints()-1):
                acc = self.get_acc(robot, parallel_number, point_number)
                if (abs(acc) > robot.GetDOFMaxAccel()).any():
                    acc_not_valid.append([parallel_number, point_number])
        return acc_not_valid

    def get_failed_vel(self, robot):
        """ Method returns a list [parallel_number, point_number] (Path.data structure) with all velocities above
        the limit.
        Args:
            robot: (Robot) is the robot object.
        Returns:
            List int[m][2] [parallel_number, point_number].
        Examples:
            >>> vel_fail = path.get_failed_vel(robot)
        """

        vel_not_valid = []
        for parallel_number in range(len(self.data)):
            for point_number in range(1,self.data[parallel_number].GetNumWaypoints()-1):
                vel = self.get_velocity(robot, parallel_number, point_number)
                if (vel > robot.GetDOFMaxVel()).any():
                    vel_not_valid.append([parallel_number, point_number])
        return vel_not_valid

    def plot_velocities(self, turbine, parallel_number):
        """ Plot joint velocities of specific parallel.
        Args:
            turbine: (@ref Turbine) is the turbine object.
            parallel_number: (int) is the parallel number to get the joint.
        Returns:
            velocity graphic.
        Examples:
            >>> path.plot_velocities(turbine, 0)
        """

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
        """ Plot joint accelerations of specific parallel.
        Args:
            turbine: (@ref Turbine) is the turbine object.
            parallel_number: (int) is the parallel number to get the joint.
        Returns:
            acceleration graphic.
        Examples:
            >>> path.plot_acc(turbine, 0)
        """

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

    def plot_pos(self, turbine, parallel_number):
        """ Plot joint position of specific parallel.

        Args:
            turbine: (@ref Turbine) is the turbine object.
            parallel_number: (int) is the parallel number to get the joint.

        Returns:
            position graphic.

        Examples:
            >>> path.plot_pos(turbine, 0)
        """

        pos = []
        dtimes = []
        min_pos = turbine.robot.GetActiveDOFLimits()[0]
        max_pos = turbine.robot.GetActiveDOFLimits()[1]

        N = self.data[parallel_number].GetNumWaypoints()

        for i in range(N):
            pos.append(self.get_joint(turbine.robot, parallel_number, i))
            dtimes.append(self.get_deltatime(parallel_number, i))

        pos = array(pos)
        f, ax = plt.subplots(turbine.robot.GetDOF(), sharex=True)
        dtimes = cumsum(dtimes)

        for i in range(turbine.robot.GetDOF()):
            ax[i].plot(dtimes, pos[:, i], color='b')
            ax[i].plot(dtimes, ones(N) * max_pos[i], color='r')
            ax[i].plot(dtimes, ones(N) * min_pos[i], color='r')
            ax[i].set_title('Joint Position ' + str(i))

        f.subplots_adjust(hspace=0.3)
        plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
        plt.show()
        return

    def plot_compare_mls_dijkstra(self, turbine, parallel_number):
        """ Plot joint position of specific parallel.

        Args:
            turbine: (@ref Turbine) is the turbine object.
            parallel_number: (int) is the parallel number to get the joint.

        Returns:
            position graphic.

        Examples:
            >>> path.plot_compare_mls_dijkstra(turbine, 0)
        """

        pos = []
        dtimes = []
        min_pos = turbine.robot.GetActiveDOFLimits()[0]
        max_pos = turbine.robot.GetActiveDOFLimits()[1]

        N = self.data[parallel_number].GetNumWaypoints()

        for i in range(N):
            pos.append(self.get_joint(turbine.robot, parallel_number, i))
            dtimes.append(self.get_deltatime(parallel_number, i))

        pos = array(pos)
        f, ax = plt.subplots(turbine.robot.GetDOF(), sharex=True)
        dtimes = cumsum(dtimes)

        for i in range(turbine.robot.GetDOF()):
            ax[i].plot(dtimes, pos[:, i], color='b')
            ax[i].plot(dtimes, ones(N) * max_pos[i], color='r')
            ax[i].plot(dtimes, ones(N) * min_pos[i], color='r')
            ax[i].plot(dtimes, array(self.joint_dijkstra[parallel_number])[:, i], 'go')
            ax[i].set_title('Joint Position ' + str(i))

        f.subplots_adjust(hspace=0.3)
        plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
        plt.show()
        return

    def plot_failed_points(self, robot, vis):
        manip = robot.GetActiveManipulator()
        valid_points = []
        vel_not_valid = []
        acc_not_valid = []

        with robot:
            for parallel_number in range(len(self.data)):
                N = self.data[parallel_number].GetNumWaypoints()
                for point_number in range(N):
                    joint = self.get_joint(robot, parallel_number, point_number)
                    vel = self.get_velocity(robot, parallel_number, point_number)
                    acc = self.get_acc(robot, parallel_number, point_number)
                    robot.SetDOFValues(joint)
                    pos = manip.GetTransform()[0:3,3]

                    if point_number > 0 and point_number < N-1:
                        if (vel > robot.GetDOFMaxVel()).any():
                            vel_not_valid.append(pos)
                            continue

                        if (acc > robot.GetDOFMaxAccel()).any():
                            acc_not_valid.append(pos)

                    valid_points.append(pos)

        _  = vis.plot(valid_points, 'valid', (0, 0, 0))
        _ = vis.plot(vel_not_valid, 'velnot', (1, 0, 0))
        _ = vis.plot(acc_not_valid, 'accnot', (0, 1, 0))
        return

    def compute_iksol(self, turbine, organized_rays_list, interpolation):
        """ Given organized_rays_list (cartesian points x-y-z-nx-ny-nz), this method returns robot's joints solution,
        using inverse kinematics (openrave ikfast).

        Args:
            turbine: (@ref Turbine) turbine object
            organized_rays_list: (float[m][n<SUB>i</SUB>][6]) zigzagging parallels
            interpolation: (float) distance between points in parallels
        Returns:
             joint values (all possibilitites)
        Examples:
            >>> iksols, dtimes = path.move_dijkstra(turbine, rays_list, 3e-3)
        """

        self.organized_rays_list = mathtools.equally_spacer(organized_rays_list, interpolation)
        self.iksols = []

        for i, organized_rays in enumerate(organized_rays_list):
            iksol = robot_utils.compute_robot_joints(turbine, organized_rays, True)
            self.iksols.append(iksol)
        return

    @staticmethod
    def move_dijkstra(turbine, iksols):
        """ Given iksols (joint solutions) and dtimes, this method returns robot's joints solution, using
        the Dijktra algorithm for path optimization: @see https://en.wikipedia.org/wiki/Dijkstra%27s_algorithm

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
        dtimes_list = []

        for i, iksol in enumerate(iksols):

            iksol = copy(iksol)

            # Add the first element of the next parallel to the end of the current parallel
            # if the next parallel exists (current is not the last)
            if i != len(iksols) - 1:
                iksol.append(iksols[i + 1][0])

            # Add the ik solutions of the previous parallel to the beginning of the current parallel
            # if the previous parallel exists (current is not the first)
            if i != 0:
                iksol.insert(0, [joint_path_list[-1][-1]])

            # Compute deltatimes from operational space
            dtimes = compute_dtimes_from_joints(turbine, [j[0] for j in iksol])

            # Remove the last joint from planning
            iksol = array([array(j)[:, :-1] for j in iksol])

            # Call dijkstra
            joint_path, path, min_cost, adj, cost = dijkstra_acc.make_dijkstra(iksol, dtimes, vel_limits, acc_limits,
                                                                               True)

            # Remove the last solution of the current parallel because it is the first solution of the next parallel
            # if the current parallel is not the last one
            if i != len(iksols) - 1:
                joint_path = joint_path[:-1]
                dtimes = dtimes[:-1]

            # Remove the first solution of the current parallel because it is the last solution of the previous parallel
            # if the current parallel is not the first one
            if i != 0:
                joint_path = joint_path[1:]
                dtimes = dtimes[1:]

            dtimes_list.append(dtimes)

            # Add the last joint again, with zero value
            joint_path_complete = zeros((len(joint_path), robot.GetDOF()))
            joint_path_complete[:len(joint_path), :len(joint_path[0])] = joint_path
            joint_path_list.append(joint_path_complete)

        return joint_path_list, dtimes_list

    def break_grid(self, turbine):
        init = 0
        grid_break = []
        robot = turbine.robot

        while True:
            acc_fail = self.get_failed_acc(robot)
            vel_fail = self.get_failed_vel(robot)
            joint_dijkstra = []

            if len(vel_fail) > 0 or len(acc_fail) > 0:
                if len(vel_fail) > 0:
                    first_vel_fail_parallel = vel_fail[0][0]
                else:
                    first_vel_fail_parallel = len(self.data)
                if len(acc_fail) > 0:
                    first_acc_fail_parallel = acc_fail[0][0]
                else:
                    first_acc_fail_parallel = len(self.data)
                fail_parallel = min(first_vel_fail_parallel, first_acc_fail_parallel)
                grid_break.append(fail_parallel)

                if init == fail_parallel:
                    return grid_break, False

                if len(grid_break) > 3:
                    return grid_break, False

                for i in range(0,fail_parallel):
                    joint_dijkstra.append(self.joint_dijkstra[i])

                joint_dijkstra2, _ = self.move_dijkstra(turbine, self.iksols[fail_parallel:])

                for i in range(len(joint_dijkstra2)):
                    joint_dijkstra.append(joint_dijkstra2[i])

                self.joint_dijkstra = joint_dijkstra

                joint, joint_velocity, joint_acc = self.mls_parallels(turbine, self.joint_dijkstra, self.dtimes)

                self.create_trajectory_data(turbine,joint,joint_velocity,joint_acc)

                init = fail_parallel
            else:
                return grid_break, True
        return grid_break, True


    def mls_parallels(self, turbine, joints_parallels, dtimes_parallels):
        """ This method calls mls_joints to smooth the trajectories (all parallels)
        Args:
            turbine: (@ref Turbine) turbine object
            joint_path: (float[m][n<SUB>i</SUB>][nDOF]) list of joints for all parallels
            dtimes_parallels: (float[m][n<SUB>i</SUB>][1]) list of dtimes for all parallels
        Returns:
            Four lists: smooth joint values, joint velocities,  joint accelerations, and deltatimes
        Examples:
            >>> joint_path, joint_velocity_path, joint_acc_path, dtimes = path.mls_parallels(turbine, joint_path, dtimes_path)
        """

        new_joints_parallels = []
        new_joints_vel_parallels = []
        new_joints_acc_parallels = []
        for joints, dtimes in zip(joints_parallels, dtimes_parallels):
            dtimes = list(dtimes)
            dt1 = dtimes[1]
            dt12 = dtimes[2] + dt1
            v = vander([dt1, dt12], 3)[:2, :2]

            j1 = joints[1]
            j0 = joints[0]
            j2 = joints[2]

            x = linalg.solve(v, [j1 - j0, j2 - j0])
            nv = [1, -1] * v
            jm1, jm2 = dot(nv, x) + j0


            dtn1 = dtimes[-1]
            dtn12 = dtimes[-2] + dt1
            v = vander([-dtn1, -dtn12], 3)[:2, :2]

            jnm1 = joints[-2]
            jn = joints[-1]
            jnm2 = joints[-3]

            x = linalg.solve(v, [jnm1 - jn, jnm2 - jn])
            nv = [1, -1] * v
            jnm1, jnm2 = dot(nv, x) + jn


            joints = array([jm2, jm1] + list(joints) + [jnm1, jnm2])
            dtimes = array([dtimes[0], dtimes[2], dtimes[1]] + dtimes[1:] + [dtimes[-1],dtimes[-2]])
            new_joints, new_vel, new_acc = self.mls_joints(
                turbine, joints, error=5.e-3, mls_degree=6, times=cumsum(dtimes))
            new_joints_parallels += [new_joints[2:-2]]
            new_joints_vel_parallels += [new_vel[2:-2]]
            new_joints_acc_parallels += [new_acc[2:-2]]

        return new_joints_parallels, new_joints_vel_parallels, new_joints_acc_parallels

    def mls_joints(self, turbine, joints, error, mls_degree = 6, times=None):
        """ The Dijkstra optimization method will calculate the best path for the given discretization. Discretization
        will generate some "jumps", and infeasible velocities/accelerations. A moving least square method was developed
        for path smoothness. The MLS generates a trajectory that will not pass by all waypoints, but it will choose the
        best scale given the maximum required error in cartesian space (meters).
        It returns joints, joint velocities, accelerations.
        This method chooses the best scale given the maximum required error in cartesian space.
        It receives one parallel.
        Args:
            turbine: (@ref Turbine) turbine object
            joints: (float[n<SUB>i</SUB>][nDOF]) list of joints]
            error: (float) maximum required error in cartesian space (meters)
            scale: (float) MLS scale
            times: (float[n<SUB>i</SUB>]) cumsum(deltatimes)
        Returns:
            lists: smooth joint values, joint velocities,  joint accelerations
        Examples:
            >>> joints, joint_velocity, joint_acc, dtimes = path.mls_joints(turbine, joints)
        """

        if times is None:
            times = cumsum(compute_dtimes_from_joints(turbine, joints))

        scale = (((times[2:] - times[1:-1]) ** 2 + (times[1:-1] - times[:-2]) ** 2) / 2) ** 0.5
        scale = array([abs(times[1]-times[0])] + list(scale) + [abs(times[-1]-times[-2])])

        counter = 0
        while True:
            counter += 1
            new_joints = []
            new_joints_velocities = []
            new_joints_acc = []
            for i, joint in enumerate(array(joints).T):
                j, v, a = mathtools.legMLS(joint, times, times, mls_degree, scale)
                new_joints += [j]
                new_joints_velocities += [v]
                new_joints_acc += [a]
            new_joints = array(new_joints).T
            new_joints_velocities = array(new_joints_velocities).T
            new_joints_acc = array(new_joints_acc).T
            ans, actual_error = self.joint_error(turbine.robot, joints, new_joints, error)
            if ans:
                if actual_error < error * 1e-2:
                    scale *= 2
                    if counter > 50:
                        break
                    continue
                else:
                    break
            scale -= .05 * scale
        return array(new_joints), array(new_joints_velocities), array(new_joints_acc)

    @staticmethod
    def joint_error(robot, joints_a, joints_b, error_max):
        """ This method calculates the error between two lists of joints (a and b) in cartesian
        space. It iterates in lists joints_b and joints_a, setting the robot DOF values and comparing the end-effector
        position. Orientations are not compared.
        Args:
            robot: (Robot) object
            joints_a: (float[n<SUB>i</SUB>][nDOF]) list of joints_a
            joints_b: (float[n<SUB>i</SUB>][nDOF]) list of joints_b to be compared with.
        Returns:
            mean(error) + 0.5*std (float) and new_points
        Examples:
            >>> e,_ = path.joint_error(robot, joints_a, joints_b)
        """

        rays_a = get_ee_ray(robot, joints_a)
        pos_b = get_ee_positions(robot, joints_b)
        errors = []
        for j in range(len(rays_a)):
            x1 = array(rays_a[j][0:3])
            x2 = x1 + array(rays_a[j][3:6])
            error = mathtools.distance_point_line_3d(x1, x2, array(pos_b[j]))
            if error > error_max:
                return False, error
            errors.append(error)
        return True, max(errors)


    @staticmethod
    def parallels_transitions(turbine, joints_parallel, joints_vel_parallel, joints_acc_parallel, times_parallel,
                              step=.1, number_of_points=8.) :
        """ Given all parallels solutions, i.e., joints solutions for all waypoints split in parallels (lists),
        this method computes cubic polynomials to interpolate the two points: end of a parallel - begin of the next
        parallel, in joint space, considering velocities. It will return the parallels and the computed transition
        paths between them.
        Obs.: This is taking too long. No improvements are being made since plates can be used on the transition part.
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
            vel_0 = joints_vel_parallel[2*i][-1]
            joint_1 = joints_parallel[2*i + 1][0]
            vel_1 = joints_vel_parallel[2*i + 1][0]

            c = mathtools.legn_path(3, [joint_0, joint_1], [vel_0, vel_1])
            joints = []
            joints_vel = []
            joints_acc = []
            times = []
            dt = min([step, 1. / number_of_points])
            for t in linspace(0, 1., max([1. / step, number_of_points])):
                joints.append(legendre.legval(t, c))
                joints_vel.append(legendre.legval(t, legendre.legder(c, 1)))
                joints_acc.append(legendre.legval(t, legendre.legder(c, 2)))
                times.append(dt)

            joints = joints[1:-1]
            joints_vel = joints_vel[1:-1]
            joints_acc = joints_acc[1:-1]
            times_parallel[2 * i + 1][0] = times[-1]
            times = times[1:-1]

            joints_parallel.insert(2 * i + 1, joints)
            joints_vel_parallel.insert(2 * i + 1, joints_vel)
            joints_acc_parallel.insert(2 * i + 1, joints_acc)
            times_parallel.insert(2 * i + 1, times)

        return joints_parallel, joints_vel_parallel, joints_acc_parallel, times_parallel

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

def get_ee_positions(robot, joints):
    ee_pos = []
    manip = robot.GetActiveManipulator()
    with robot:
        for joint in joints:
            robot.SetDOFValues(joint)
            ee_pos.append(manip.GetTransform()[0:3,3])
    return ee_pos

def get_ee_orientation(robot, joints):
    ee_ori = []
    manip = robot.GetActiveManipulator()
    with robot:
        for joint in joints:
            robot.SetDOFValues(joint)
            ee_ori.append(manip.GetTransform()[0:3,0:3])
    return ee_ori

def get_ee_ray(robot, joints):
    ee_ray = []
    manip = robot.GetActiveManipulator()
    with robot:
        for joint in joints:
            robot.SetDOFValues(joint)
            t = manip.GetTransform()
            ee_ray.append([t[0,3], t[1,3], t[2,3], t[0,0], t[1,0], t[2,0]])
    return ee_ray
