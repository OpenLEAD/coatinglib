from numpy import array, linalg, dot, zeros, inf, vstack, mean, std, cumsum, abs, sign, ones
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

    def execute(self, turbine, threshold=5e-2, dtimes = None):
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

        joint_path, rays_list = self.move_dijkstra(turbine, self.rays, threshold)
        if len(joint_path) == 0:
            return
        # joint_path = self.refine_dijkstra(turbine, joint_path, rays_list, threshold)

        new_joint_path, new_joint_velocity_path, new_joint_acc_path, new_times_path = self.smooth_joint_MLS(
            turbine, joint_path, dtimes)

        # new_joint_path, new_joint_velocity_path, new_joint_acc_path, new_times_path = \
        #    self.replanning(turbine, new_joint_path, new_joint_velocity_path, new_joint_acc_path, new_times_path)

        acc = abs(vstack(new_joint_acc_path))
        for i,vel_paralell in enumerate(new_joint_velocity_path):
            for j,vel in enumerate(vel_paralell):
                if (vel[:-1] > turbine.robot.GetDOFMaxVel()[:-1]).any():
                    print '(',i, ',', j,')', '| parallel length = ', str(j),'/',len(vel_paralell)-1, '| parallel number = ', str(i),'/',len(new_joint_velocity_path)-1
                    print 'vel =', vel


        if (acc > turbine.robot.GetDOFMaxAccel()).any():
            print 'acc max fail'

        for i,acc_paralell in enumerate(new_joint_acc_path):
            for j,acc in enumerate(acc_paralell):
                if (acc[:-1] > turbine.robot.GetDOFMaxVel()[:-1]).any():
                    print '(',i, ',', j,')', '| parallel length = ', str(j),'/',len(acc_paralell)-1, '| parallel number = ', str(i),'/',len(new_joint_acc_path)-1
                    print 'acc =', acc

        self.success = True

        ind = str()
        for i in range(turbine.robot.GetDOF()): ind += str(i) + ' '

        cs = ConfigurationSpecification()
        _ = cs.AddGroup('joint_values ' + turbine.robot.GetName() + ' ' + ind, len(turbine.robot.GetActiveDOFIndices()),
                        'cubic')
        cs.AddDerivativeGroups(1, False)
        cs.AddDerivativeGroups(2, False)
        _ = cs.AddDeltaTimeGroup()

        TRAJ = []
        for i in range(len(new_joint_path)):
            traj = RaveCreateTrajectory(turbine.env, '')
            traj.Init(cs)
            for j in range(len(new_joint_path[i])):
                waypoint = list(new_joint_path[i][j])
                waypoint.extend(list(new_joint_velocity_path[i][j]))
                waypoint.extend(list(new_joint_acc_path[i][j]))
                waypoint.extend([new_times_path[i][j]])
                traj.Insert(traj.GetNumWaypoints(), waypoint)
            TRAJ.append(traj)
        self.data = TRAJ
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
            ax[i].plot(dtimes,velocities[:,i], color='b')
            ax[i].plot(dtimes, ones(N)*max_vel[i], color='r')
            ax[i].plot(dtimes, -ones(N)*max_vel[i], color = 'r')
            ax[i].set_title('Joint Velocity '+str(i))

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
            ax[i].plot(dtimes,acc[:,i], color='b')
            ax[i].plot(dtimes, ones(N)*max_acc[i], color='r')
            ax[i].plot(dtimes, -ones(N)*max_acc[i], color = 'r')
            ax[i].set_title('Joint Acc '+str(i))

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
        organized_rays_list = mathtools.equally_spacer(organized_rays_list, interpolation)
        true_deep = False
        t = Time.time()
        for i,organized_rays in enumerate(organized_rays_list):

            for deep in [True]:
                if deep == False and true_deep == True:
                    continue
                try:
                    joints = planning.joint_planning(turbine, organized_rays, deep)
                    total_joints = []
                    for i in range(len(joints)-1):
                        total_joints.append(len(joints[i])*len(joints[i+1]))
                    V = sum(total_joints)+len(joints[0])+len(joints[-1])
                    print 'number of vertices for dijkstra = ', V
                except IndexError:
                    true_deep = True
                    continue
                if len(joint_path_list) != 0:
                    joints.insert(0, [joint_path_list[-1][-1]])
                    dtimes = compute_dtimes_from_joints(turbine, [j[0] for j in joints])
                    joints = array([array(j)[:,:-1] for j in joints])
                    joint_path, path, min_cost, adj, cost = planning.make_dijkstra(joints, dtimes, vel_limits, acc_limits, True)

                    print min_cost
                    if min_cost != inf:
                        joint_path_complete = zeros((len(joint_path),robot.GetDOF()))
                        joint_path_complete[:len(joint_path), :len(joint_path[0])] = joint_path
                        joint_path_list.append(joint_path_complete[1:])
                        break

                else:
                    dtimes = compute_dtimes_from_rays(turbine, organized_rays)
                    joints = array([array(j)[:,:-1] for j in joints])
                    joint_path, path, min_cost, adj, cost = planning.make_dijkstra(joints, dtimes, vel_limits, acc_limits, True)
                    print  min_cost
                    if min_cost != inf:
                        joint_path_complete = zeros((len(joint_path),robot.GetDOF()))
                        joint_path_complete[:len(joint_path), :len(joint_path[0])] = joint_path
                        joint_path_list.append(joint_path_complete[1:])
                        break
            else:
                return [], organized_rays_list
        print 'time = ', Time.time() - t
        return joint_path_list, organized_rays_list

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

    def smooth_joint_MLS(self, turbine, joint_path, dtimes=None):
        """ The Dijkstra optimization method will calculate the best path for the given discretization. Discretization
        will generate some "jumps", and infeasible velocities/accelerations. A moving least square method was developed
        for path smoothness.

        Args:
            turbine: (@ref Turbine) turbine object
            joint_path: (float[m][n<SUB>i</SUB>][nDOF]) list of joints for all parallels

        Returns:
            Four lists: smooth joint values, joint velocities,  joint accelerations, and deltatimes

        Examples:
            >>> joint_path, joint_velocity_path, joint_acc_path, dtimes = path.smooth_joint_MLS(turbine, joint_path)
        """

        new_joint_path = []
        new_joint_velocity_path = []
        new_joint_acc_path = []
        new_dtimes_path = []

        for joints in joint_path:
            scale = 3.
            try:
                self.joints_MLS(turbine, joints, scale, dtimes)
            except TypeError:
                print 'jnts shape3 - ', array(joints).shape
                print 'scale 3- ', scale


            try:
                new_joints, new_joints_velocities, new_joints_acc, new_dtimes = self.joints_MLS(turbine, joints, scale, dtimes)
            except TypeError:
                print 'jnts shape2 - ', array(joints).shape
                print 'scale 2- ', scale
                raise

            new_joint_path += [new_joints]
            new_dtimes_path += [new_dtimes]
            new_joint_velocity_path += [new_joints_velocities]
            new_joint_acc_path += [new_joints_acc]

        return new_joint_path, new_joint_velocity_path, new_joint_acc_path, new_dtimes_path

    def joints_MLS(self, turbine, joints, scale=3., times=None):

        if times is None:
            times = cumsum(compute_dtimes_from_joints(turbine, joints))

        while scale > 0:
            new_joints = []
            new_joints_velocities = []
            new_joints_acc = []
            for i, joint in enumerate(array(joints).T):
                j, v, a = mathtools.legMLS(joint, times, times, 6, scale)
                new_joints += [j]
                new_joints_velocities += [v]
                new_joints_acc += [a]
            new_joints = array(new_joints).T
            new_dtimes = compute_dtimes_from_joints(turbine, new_joints)
            new_joints_velocities = array(new_joints_velocities).T
            new_joints_acc = array(new_joints_acc).T
            error, _ = self.joint_error(turbine.robot, joints, new_joints)
            if error <= 2.5e-3:  # HARDCODED
                break
            scale -= .1
        return array(new_joints), array(new_joints_velocities), array(new_joints_acc), array(new_dtimes)

    @staticmethod
    def joint_error(robot, joints, new_joints):
        manip = robot.GetActiveManipulator()
        error = []
        new_points = []
        with robot:
            for j in range(len(new_joints)):
                robot.SetDOFValues(new_joints[j])
                P0 = manip.GetTransform()[0:3, 3]
                robot.SetDOFValues(joints[j])
                P1 = manip.GetTransform()[0:3, 3]
                error.append(linalg.norm(P0 - P1))
                new_points += [P0]
        return mean(error) + .5 * std(error), new_points

    def replanning(self, turbine, new_joint_path, new_joint_velocity_path, new_joint_acc_path, new_dtimes_path):

        for k in range(len(new_joint_path)):
            q, dq, ddq, dt = array(new_joint_path[k]), array(new_joint_velocity_path[k]), array(new_joint_acc_path[k]), array(new_dtimes_path[k])

            for j in range(len(ddq[0])):
                timeshift_list = []
                # start_h = -1
                # start = -1
                max_accel = turbine.robot.GetDOFMaxAccel()[j]
                for i, ddqij in enumerate(ddq[:,j]):
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

                    if abs(ddqij) > max_accel*0.95:
                        timeshift_list.append(i)

                if len(timeshift_list) > 0:
                    q[:,j], dt = self.timeshift(timeshift_list, max_accel*0.95, q[:,j], dq[:,j], ddq[:,j], dt)
                    try:
                        q, dq, ddq, dt = self.joints_MLS(turbine, q, 3., cumsum(dt))
                    except TypeError:
                        print 'dt type - ', type(dt)

            new_joint_path[k], new_joint_velocity_path[k], new_joint_acc_path[k], new_dtimes_path[k] = q, dq, ddq, dt
        return new_joint_path, new_joint_velocity_path, new_joint_acc_path, new_dtimes_path

    @staticmethod
    def timeshift(timeshift_list, max_accel, q, dq, ddq, dt):
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
            if t == len(q)-1:
                continue

            dq0 = dq[t]
            delta_q = q[t+1] - q[t]
            limit = sign(ddq[t]) * max_accel * 0.95
            tf = (-dq0 + (dq0 ** 2 + 2 * delta_q * limit) ** 0.5) / limit

            print 't - ', t
            print 'ddq[t] - ', ddq[t]
            print 'limit - ', limit
            print 'dq0 - ', dq0
            print 'dqf - ', dq[t+1]
            print 'delta_q - ', delta_q
            print 'tf - ', tf
            print 'dt[t+1] -', dt[t + 1]

            dt[t+1] = tf

        return q, dt

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
    """ Given the real blade angle:
    1) rotate the blades (update the environment);
    2) organize trajectories (removing empty points, adding borders,
    and making zigzagging lists);
    3) Dijkstra algorithm.

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
    """ Given the joints solution, this method access the config file (velocity requirements) and compute the required
    delta times between joints.

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

    return compute_dtimes_from_rays(turbine,rays)


def compute_dtimes_from_rays(turbine,rays):

    rays = array(rays)
    v = turbine.config.coating.coating_speed

    return array([0.]+list(linalg.norm(rays[1:,0:3]-rays[:-1,0:3],axis=1)/v))


def jusante_grids():
    """ Return the int numbers from jusante grids.

    Returns:
        grid numbers

    Examples:
        >>> grid_nums = jusante_grids()
    """

    grid_nums = range(0, 15)
    grid_nums.extend(range(17, 20))
    grid_nums.extend(range(22, 25))
    grid_nums.extend(range(67, 70))
    grid_nums.extend(range(72, 77))
    grid_nums.extend(range(78, 80))
    return grid_nums


def montante_grids():
    """ Return the int numbers from montante grids.

    Returns:
        grid numbers

    Examples:
        >>> grid_nums = montante_grids()
    """

    grid_nums = range(30, 50)
    grid_nums.extend(range(51, 55))
    grid_nums.extend(range(56, 60))
    grid_nums.extend([77])
    return grid_nums


def lip_grids():
    """ Return the int numbers from lip grids.

    Returns:
        grid numbers

    Examples:
        >>> grid_nums = lip_grids()
    """

    return [0, 1, 2]


def border_grids():
    """ Return the int numbers from border grids.

    Returns:
        grid numbers

    Examples:
        >>> grid_nums = border_grids()
    """

    grid_nums = range(60, 65)
    grid_nums.extend(range(25, 29))
    grid_nums.extend([85])
    return grid_nums