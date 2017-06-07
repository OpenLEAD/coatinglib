from numpy import array, linalg, dot, zeros, inf, vstack, mean, std, cumsum
from numpy import max as pymax
from numpy import sum as pysum
import db
import planning
from openravepy import ConfigurationSpecification, interfaces, planningutils, RaveCreateTrajectory
import mathtools
import time
import errno
from scipy.linalg import logm, expm
from math import acos
import rail_place
from xml.etree import ElementTree as ET
from os import listdir, makedirs
from os.path import realpath, splitext, join

class path:
    """
    Class to compute, store, and serialize joint path.
    Robot's full trajectories are stored in OpenRave format Trajectory Class, as below:
    
    joint_values
    joint_vel_calues
    joint_acc_values
    deltatime
    """
    
    def __init__(self, rays=None):
        """
        Constructor fo class path.

        @param rays (nx6 np.array) are the points and normal vectors of the points to coat (x,y,z,nx,ny,nz).
        """
        self.rays = rays
        self.data = []
        self.success = False

    def execute(self, turbine, blade, threshold=5e-2):
        """
        Method to compute joint_values, joint_velocities, joint_accelerations and deltatimes.
        It uses the Dijkstra planning algorithm (see move_dijkstra function).
        
        @param turbine is the turbine object.
        @param blade is the blade object (you may use DB.load_blade(), for instance)
        @param threshold (float) is the interpolation threshold, as rays are usually well spaced.
        """
        
        if self.rays is None: return
        
        joint_path, rays_list = move_dijkstra(turbine, blade, self.rays, threshold)
        if len(joint_path) == 0:
            return
        joint_path = refine_dijkstra(turbine, joint_path, rays_list, threshold)
        
        new_joint_path, new_joint_velocity_path, new_joint_acc_path, new_times_path = smooth_joint_MLS(turbine, joint_path)
        maxAccel = turbine.robot.GetDOFMaxAccel()
        maxVel = turbine.robot.GetDOFMaxVel()
        vel = vstack(new_joint_velocity_path)
        acc = vstack(new_joint_acc_path)
        if (vel>maxVel).any(): return
        if (acc>maxAccel).any(): return
        self.success = True

        ind = str()
        for i in range(turbine.robot.GetDOF()): ind+=str(i)+' '

        cs = ConfigurationSpecification()
        _ = cs.AddGroup('joint_values '+turbine.robot.GetName()+' '+ind, len(turbine.robot.GetActiveDOFIndices()), 'cubic')
        cs.AddDerivativeGroups(1,False)
        cs.AddDerivativeGroups(2,False)
        _ = cs.AddDeltaTimeGroup()

        TRAJ = []
        for i in range(len(new_joint_path)):
            traj = RaveCreateTrajectory(turbine.env,'')
            traj.Init(cs)
            for j in range(len(new_joint_path[i])):
                waypoint = list(new_joint_path[i][j])
                waypoint.extend(list(new_joint_velocity_path[i][j]))
                waypoint.extend(list(new_joint_acc_path[i][j]))
                waypoint.extend([new_times_path[i][j]])
                traj.Insert(traj.GetNumWaypoints(),waypoint)
            TRAJ.append(traj)
        self.data = TRAJ
        return

    def serialize(self, directory=''):
        """
        Method to serialize data in OpenRave format output is an xml file.
        @param directory (str) is the relative path to the folder.
        """
        
        try:
            makedirs(directory)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
            
        for i, traj in enumerate(self.data):
            tree = ET.XML(traj.serialize())
            with open(join(directory,'trajectory_'+str(i)+'.xml'), 'w') as f:
                f.write(ET.tostring(tree))
        return

    def deserialize(self, turbine, directory):
        """
        Method to deserialize data in OpenRave format.
        @param turbine is the turbine object.
        @param directory (str) is the relative path to the folder.
        """
        
        cs = ConfigurationSpecification()
        _ = cs.AddGroup('joint_values', len(turbine.robot.GetActiveDOFIndices()), 'linear')
        cs.AddDerivativeGroups(1,False)
        cs.AddDerivativeGroups(2,False)
        _ = cs.AddDeltaTimeGroup()
        
        path = os.realpath(path)
        TRAJ = []
        onlyfiles = [f for f in listdir(directory) if isfile(join(directory, f))]
        for afile in onlyfiles:
            filename, file_extension = splitext(afile)
            if file_extension == '.xml':
                traj = RaveCreateTrajectory(turbine.env,'')
                traj.Init(cs)
                xml = ET.parse(open(join(directory,afile)))
                traj.deserialize(ET.tostring(xml.getroot()))
            TRAJ.append(traj)
        self.data = TRAJ
        return

    def simulate(self, robot, parallel_number):
        """
        Method to simulate and visualize, in real time, the robot performing the coat.
        @param robot is the robot object.
        @param parallel_number (int) is parallel to be simulated.
        """
        
        ans = robot.GetController().SetPath(path.data[parallel_number])

    def get_joint(self, robot, parallel_number, point_number):
        """
        Method to get a specific joint_value from path.
        @param robot is the robot object.
        @param parallel_number (int) is the parallel number to get the joint.
        @param point_number (int) is the specific point in the parallel to get the joint _value.
        """
        traj = path.data[parallel_number]
        spec = traj.GetConfigurationSpecification()
        return spec.ExtractJointValues(traj.GetWaypoint(point_number),robot,range(robot.GetDOF()))

    def get_velocity(self, robot, parallel_number, point_number):
        """
        Method to get a specific joint_velocity from path.
        @param robot is the robot object.
        @param parallel_number (int) is the parallel number to get the joint_velocity.
        @param point_number (int) is the specific point in the parallel to get the joint_velocity.
        """
        traj = path.data[parallel_number]
        spec = traj.GetConfigurationSpecification()
        return spec.ExtractJointValues(traj.GetWaypoint(point_number),robot,range(robot.GetDOF()),1)

    def get_acc(self, robot, parallel_number, point_number):
        """
        Method to get a specific joint_acceleration from path.
        @param robot is the robot object.
        @param parallel_number (int) is the parallel number to get the joint_acceleration.
        @param point_number (int) is the specific point in the parallel to get the joint_acceleration.
        """
        traj = path.data[parallel_number]
        spec = traj.GetConfigurationSpecification()
        return spec.ExtractJointValues(traj.GetWaypoint(point_number),robot,range(robot.GetDOF()),2)

    def get_deltatime(self, robot, parallel_number, point_number):
        """
        Method to get a specific deltatime from path.
        @param robot is the robot object.
        @param parallel_number (int) is the parallel number to get the deltatime.
        @param point_number (int) is the specific point in the parallel to get the deltatime.
        """
        traj = path.data[parallel_number]
        spec = traj.GetConfigurationSpecification()
        return spec.ExtractDeltaTime(traj.GetWaypoint(point_number))
                                       
        

def organize_rays_in_parallels(DB, grid):
    rays = DB.compute_rays_from_grid(grid)

    organized_rays = []
    for i in range(0,len(rays)):
        not_empty_rays = mathtools.notempty(rays[i])
        if len(not_empty_rays)==0:
            continue
        if i%2==0:  
            organized_rays.append(not_empty_rays)
        else:
            organized_rays.append(list(reversed(not_empty_rays)))
    return organized_rays

def base_grid_validation(turbine, psa, DB, grid, threshold = 5e-2):
    """
    Given the real blade angle:
    1) rotate the blades (update the environment);
    2) organize trajectories (removing empty points, adding borders,
    and making zigzagging lists);
    3) Dijkstra algorithm.

    @param turbine is the turbine object.
    @param psa (tuple 1x3) primary, secondary, alpha (base position)
    @param DB database object.
    @param grid (int) grid to be coated.
    @param threshold (float) is the interpolation threshold, as rays are usually well spaced.
    """
    
    robot = turbine.robot
    manip = robot.GetActiveManipulator()
    
    T = DB.T
    T = dot(T,linalg.inv(turbine.blades[3].GetTransform())) #HARDCODED

    for blade in turbine.blades:
        T0 = blade.GetTransform()
        blade.SetTransform(dot(T,T0))
     
    rp = rail_place.RailPlace(psa)
    turbine.place_rail(rp)
    turbine.place_robot(rp)
    organized_rays_list = organize_rays_in_parallels(DB, grid)
    trajectory = path(organized_rays_list)
    trajectory.execute(turbine, DB.load_blade(), threshold)
    return trajectory

def generate_linear_interpolation_rays(organized_rays, blade, threshold):
    new_rays = []
    new_rays.append(organized_rays[0])
    model = blade.select_model(organized_rays[0])
    organized_rays = mathtools.filter_trajectory(organized_rays, threshold)
    for i in range(0,len(organized_rays)-1):
        points, d = mathtools.linear_interpolation_points(organized_rays[i][0:3], organized_rays[i+1][0:3], threshold)
        new_rays.extend(points[1:])
    for i in range(0,len(new_rays)):
        new_rays[i] = blade.compute_ray_from_point(new_rays[i], model)
    return new_rays

def move_dijkstra(turbine, blade, organized_rays_list, interpolation):
    robot = turbine.robot
    joint_path_list = []
    time = interpolation/turbine.config.coating.coating_speed
    limits = robot.GetDOFVelocityLimits()*time
    deep = False
    rays = []
    
    for organized_rays in organized_rays_list:
        linear_interpolation = generate_linear_interpolation_rays(organized_rays, blade, interpolation)
        rays.append(linear_interpolation)
        for deep in [False,True]:          
            try:
                joints = planning.joint_planning(turbine, linear_interpolation, deep)
            except IndexError:
                continue
            if len(joint_path_list)!=0:
                joints.insert(0,[joint_path_list[-1][-1]])
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
            return [], rays
        
    return joint_path_list, rays

def refine_dijkstra(turbine, joint_path_list, rays_list, interpolation):
    robot = turbine.robot
    time = interpolation/turbine.config.coating.coating_speed
    limits = robot.GetDOFVelocityLimits()*time
    deep = True
    new_joint_path = []
    with robot:
        for i,rays in enumerate(rays_list):
            joints_path = joint_path_list[i]
            t=0
            while True:
                new_joints = []
                t+=.01
                for j, joint in enumerate(joints_path):
                    if j==0:
                        new_joints.append([joint])
                        continue
                        
                    robot.SetDOFValues(joints_path[j-1])
                    Rx = robot.GetActiveManipulator().GetTransform()[0:3,0]
                    Rx = Rx/linalg.norm(Rx)
                    d = -dot(rays[j-1][3:6], Rx)
                    angle0 = acos(min(d,1.0))

                    robot.SetDOFValues(joint)
                    Rx = robot.GetActiveManipulator().GetTransform()[0:3,0]
                    Rx = Rx/linalg.norm(Rx)
                    d = -dot(rays[j][3:6], Rx)
                    angle1 = acos(min(d,1.0))

                    angle_tolerance_init=min([angle0,angle1])-t
                    angle_tolerance_end=max([angle0,angle1])+t
                    new_joints.append(planning.ik_angle_tolerance(turbine, rays[j],
                                                         angle_tolerance_init = angle_tolerance_init,
                                                         angle_tolerance_end = angle_tolerance_end,
                                                         number_of_phi = 24, number_of_theta = 5, deep=deep))
                if i!=0:
                    new_joints.insert(0,[new_joint_path[-1][-1]])
                    joint_path, path, min_cost, adj, cost = planning.make_dijkstra(new_joints, limits, True)
                    print min_cost
                    if min_cost != inf:
                        new_joint_path.append(joint_path[1:])
                        break
                    
                else:
                    joint_path, path, min_cost, adj, cost = planning.make_dijkstra(new_joints, limits, True)
                    print min_cost
                    if min_cost != inf:
                        new_joint_path.append(joint_path)
                        break
                    
    return new_joint_path

def smooth_trajectory(turbine, points_list, joints_list):
    new_T = []
    with turbine.robot:
        for joints in joints_list:
            T = []
            for joint in joints:
                turbine.robot.SetDOFValues(joint)
                T.append(turbine.robot.GetActiveManipulator().GetTransform())
            new_T.append(mathtools.smooth_orientation(T))
    return new_T

def smooth_joint_MLS(turbine, joint_path):
    robot = turbine.robot
    manip = robot.GetActiveManipulator()
    scale = 3.
    error = 1

    def joint_error(robot, joint_path, new_joint_path):
        error = []
        new_points = []
        with robot:
            for i in range(len(new_joint_path)):
                parallel = []
                for j in range(len(new_joint_path[i])):
                    robot.SetDOFValues(new_joint_path[i][j])
                    P0 = manip.GetTransform()[0:3,3]
                    robot.SetDOFValues(joint_path[i][j])
                    P1 = manip.GetTransform()[0:3,3]
                    error.append(linalg.norm(P0-P1))
                    parallel += [P0]
                new_points += [parallel]
        return mean(error)+.5*std(error), new_points

    while scale > 0:
        new_joint_path = []
        new_joint_velocity_path = []
        new_joint_acc_path = []
        scale-=.1
        for joints in joint_path:
                new_joints = []
                new_joints_velocities = []
                new_joints_acc = []
                for joint in array(joints).T:
                    j,v,a = mathtools.MLS(joint,array(range(len(joints))),2,scale)
                    new_joints += [j]
                    new_joints_velocities += [v]
                    new_joints_acc += [a]
                new_joint_path += [array(new_joints).T]
                new_joint_velocity_path += [array(new_joints_velocities).T]
                new_joint_acc_path += [array(new_joints_acc).T]
        error, points = joint_error(robot, joint_path, new_joint_path)
        print 'acc above 6 percent - ', pysum(vstack(new_joint_acc_path)>3,0)*6./pysum(vstack(new_joint_acc_path)>-1)
        print 'error = ', error, '| scale = ', scale, '| acc = ', pymax(vstack(new_joint_acc_path)), '| vel = ', pymax(vstack(new_joint_velocity_path))
        if error<=2.5e-3: # HARDCODED
            break
    return new_joint_path, new_joint_velocity_path, new_joint_acc_path

    

def jusante_grids():
    grid_nums = range(0,15)
    grid_nums.extend(range(17,20))
    grid_nums.extend(range(22,25))
    grid_nums.extend(range(67,70))
    grid_nums.extend(range(72,77))
    grid_nums.extend(range(78,80))
    return grid_nums

def montante_grids():
    grid_nums = range(30,50)
    grid_nums.extend(range(51,55))
    grid_nums.extend(range(56,60))
    grid_nums.extend([77])
    return grid_nums

def lip_grids():
    return [0,1,2]

def border_grids():
    grid_nums = range(60,65)
    grid_nums.extend(range(25,29))
    grid_nums.extend([85])
    return grid_nums
