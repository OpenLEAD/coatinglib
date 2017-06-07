from numpy import array, linalg, dot, zeros, inf, vstack, mean, std
from numpy import max as pymax
from numpy import sum as pysum
import db
import planning
from openravepy import ConfigurationSpecification, interfaces, planningutils
import mathtools
import time
from scipy.linalg import logm, expm
from math import acos


def organize_rays_in_parallels(DB, grid):
    DB.T = DB.turb.blades[0].GetTransform()
    parallels, borders = (DB.load_grid_to_trajectories())[grid]
    rays = DB.compute_rays_from_parallels(parallels, borders)

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

def base_grid_validation(DB, grid):
    """
    Given the real blade angle:
    1) rotate the blades (update the environment);
    2) rotate the RBF model of the blade;
    3) rotate grid points;
    4) organize trajectories (removing empty points, adding borders,
    and making one full zigzagging list);
    5) compute optimization.


def base_grid_validation_parallel(DB, grid):
    """
    Plan for each parallel.

    Keyword arguments:
    DB -- AREA DATABASE
    grid -- int
    """
    
    joint_solutions_list = []
    organized_rays = organize_rays_in_parallels(DB, grid)
    blade = DB.load_blade()
    for rays in organized_rays:
        joint_solutions_list.append(planning.compute_robot_joints(
            DB.turb, rays, 0, blade.trajectory_iter_surface))
    score = mathtools.lenlist(joint_solutions_list)*1.0/mathtools.lenlist(organized_rays)
    return score, joint_solutions_list


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
