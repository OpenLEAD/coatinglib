from numpy import array, linalg, dot, zeros, inf, vstack, mean, std, cumsum, abs, sign
import planning
from openravepy import ConfigurationSpecification, interfaces, planningutils, RaveCreateTrajectory
import mathtools
import errno
import rail_place
from xml.etree import ElementTree as ET
from os import listdir, makedirs
from os.path import realpath, splitext, join, isfile
import time as Time
from math import log

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

        vel = abs(array(new_joint_velocity_path))
        acc = abs(vstack(new_joint_acc_path))
        print 'max vel = ', turbine.robot.GetDOFMaxVel()
        for i,vel_paralell in enumerate(new_joint_velocity_path):
            for j,vel in enumerate(vel_paralell):
                if (vel[:-1] > turbine.robot.GetDOFMaxVel()[:-1]).any():
                    print '(',i, ',', j,')', '| parallel length = ', str(j),'/',len(vel_paralell)-1, '| parallel number = ', str(i),'/',len(new_joint_velocity_path)-1
                    print 'vel =', vel

        if (acc > turbine.robot.GetDOFMaxAccel()).any():
            print 'acc max fail'
            # print 'max = ', turbine.robot.GetDOFMaxAccel()
            # print acc[(acc > turbine.robot.GetDOFMaxAccel()).any(axis=1)]
            #return
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
            with open(join(directory, 'trajectory_' + str(i) + '.xml'), 'w') as f:
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

        cs = ConfigurationSpecification()
        _ = cs.AddGroup('joint_values', len(turbine.robot.GetActiveDOFIndices()), 'linear')
        cs.AddDerivativeGroups(1, False)
        cs.AddDerivativeGroups(2, False)
        _ = cs.AddDeltaTimeGroup()

        directory = realpath(directory)
        TRAJ = []
        onlyfiles = [f for f in listdir(directory) if isfile(join(directory, f))]
        for afile in onlyfiles:
            filename, file_extension = splitext(afile)
            if file_extension == '.xml':
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
            >>> for i in range(N): times.append(path.get_deltatime(turbine.robot,0,i))
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
        time = interpolation / turbine.config.coating.coating_speed
        limits = robot.GetDOFVelocityLimits() * time
        organized_rays_list = mathtools.equally_spacer(organized_rays_list, interpolation)
        true_deep = False
        for i,organized_rays in enumerate(organized_rays_list):
            for deep in [False, True]:
                if deep == False and true_deep == True:
                    continue
                try:
                    joints = planning.joint_planning(turbine, organized_rays, deep)
                except IndexError:
                    true_deep = True
                    continue
                if len(joint_path_list) != 0:
                    joints.insert(0, [joint_path_list[-1][-1]])
                    joint_path, path, min_cost, adj, cost = planning.make_dijkstra(joints, limits, True)
                    if min_cost != inf:
                        joint_path_list.append(joint_path[1:])
                        break

                else:
                    joint_path, path, min_cost, adj, cost = planning.make_dijkstra(joints, limits, True)
                    if min_cost != inf:
                        joint_path_list.append(joint_path)
                        break
            else:
                return [], organized_rays_list

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

                joint_path, path, min_cost, adj, cost = planning.make_dijkstra(new_joints, limits, True)
                if min_cost != inf:
                    new_joint_path.append(joint_path)
                else:
                    return joint_path_list
        return new_joint_path

    @staticmethod
    def smooth_joint_MLS(turbine, joint_path):
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

        robot = turbine.robot
        manip = robot.GetActiveManipulator()
        scale = 3.

        def joint_error(robot, joint_path, new_joint_path):
            error = []
            new_points = []
            with robot:
                for i in range(len(new_joint_path)):
                    parallel = []
                    for j in range(len(new_joint_path[i])):
                        robot.SetDOFValues(new_joint_path[i][j])
                        P0 = manip.GetTransform()[0:3, 3]
                        robot.SetDOFValues(joint_path[i][j])
                        P1 = manip.GetTransform()[0:3, 3]
                        error.append(linalg.norm(P0 - P1))
                        parallel += [P0]
                    new_points += [parallel]
            return mean(error) + .5 * std(error), new_points

        while scale > 0:
            new_joint_path = []
            new_joint_velocity_path = []
            new_joint_acc_path = []
            new_dtimes_path = []
            scale -= .1
            for joints in joint_path:
                new_joints = []
                new_joints_velocities = []
                new_joints_acc = []
                for joint in array(joints).T:
                    x = array(range(len(joints)))
                    j, v, a = mathtools.MLS(joint, x, x, 4, scale)
                    new_joints += [j]
                    new_joints_velocities += [v]
                    new_joints_acc += [a]

                new_joint_path += [array(new_joints).T]
                new_dtimes_path += [compute_dtimes_from_joints(turbine, new_joint_path[-1])]
                new_joint_velocity_path += [array(new_joints_velocities).T]
                new_joint_acc_path += [array(new_joints_acc).T]
            error, points = joint_error(robot, joint_path, new_joint_path)
            if error <= 2.5e-3:  # HARDCODED
                break
        return new_joint_path, new_joint_velocity_path, new_joint_acc_path, new_dtimes_path


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

    P0 = zeros(3)
    dtimes = []
    v = turbine.config.coating.coating_speed

    with turbine.robot:
        for joint in joints:
            turbine.robot.SetDOFValues(joint)
            P1 = turbine.robot.GetActiveManipulator().GetTransform()[:3, 3]
            dtimes += [linalg.norm(P1 - P0) / v]
            P0 = P1

    dtimes[0] = 0.
    return dtimes


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