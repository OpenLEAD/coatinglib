from numpy import array, linalg, dot, zeros, inf
import db
import planning
from openravepy import ConfigurationSpecification, interfaces, planningutils
import mathtools
import time
from scipy.linalg import logm, expm

def organize_rays(DB, grid):
    DB.T = DB.turb.blades[0].GetTransform()
    parallels, borders = (DB.load_grid_to_trajectories())[grid]
    rays = DB.compute_rays_from_parallels(parallels, borders)

    organized_rays = []
    for i in range(0,len(rays)):
        if i%2==0:
            organized_rays += rays[i]
        else:
            organized_rays += reversed(rays[i])
    organized_rays = mathtools.notempty(rays)
    return organized_rays

def organize_rays_in_parallels(DB, grid):
    DB.T = DB.turb.blades[0].GetTransform()
    parallels, borders = (DB.load_grid_to_trajectories())[grid]
    rays = DB.compute_rays_from_parallels(parallels, borders)

    organized_rays = []
    for i in range(0,len(rays)):
        not_empty_rays = mathtools.notempty(rays)
        if i%2==0:  
            organized_rays.append(not_empty_rays[i])
        else:
            organized_rays.append(list(reversed(not_empty_rays[i])))
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

    Keyword arguments:
    DB -- AREA DATABASE
    grid -- int
    """

    organized_rays = organize_rays(DB, grid)
    blade = DB.load_blade()
    joint_solutions = planning.compute_robot_joints_opt(DB.turb, organized_rays, 0,
                                                    blade.trajectory_iter_surface)
    score = len(joint_solutions)*1.0/len(organized_rays)

    return score, joint_solutions

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

def generate_linear_interpolation_joints(joint_solutions):
    new_joints = []
    new_joints.append(joint_solutions[0])
    for i in range(0,len(joint_solutions)-1):
        joints = mathtools.linear_interpolation_joint(joint_solutions[i], joint_solutions[i+1])
        new_joints.extend(joints[1:])
    return new_joints

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

def trajectory_generation(DB, grid):
    organized_rays = organize_rays(DB, grid)
    new_rays = generate_linear_interpolation_rays(organized_rays)
    blade = DB.load_blade()
    joint_solutions = planning.compute_robot_joints(
        DB.turb, new_rays, 0, blade.trajectory_iter_surface)
    customspec = ConfigurationSpecification()
    customspec.AddGroup('joint_values',6,'linear')
    return joint_solutions

def waitrobot(robot):
    """busy wait for robot completion"""
    while not robot.GetController().IsDone():
        time.sleep(0.01)

def movetohandposition_parallels(robot, joint_solutions_list):
    manip = robot.GetActiveManipulator()
    basemanip = interfaces.BaseManipulation(robot,plannername='birrt')
    TRAJ = []

    for joint_solutions in joint_solutions_list:
        Ts, Rs, Ps = mathtools.get_manip_transforms(robot, joint_solutions, True, True)
        robot.SetDOFValues(joint_solutions[0])
        T = []
        for i in range(0,len(Ts)-1):
            P = mathtools.linear_interpolation_points(Ps[i], Ps[i+1])
            R = mathtools.homogenous_matrix_cubic_interpolation(
                Rs[i],Rs[i+1],zeros((3,3)),zeros((3,3)),len(P))
            for j in range(0,len(P)):
                T.append(mathtools.makeTbyRP(R[j],P[j]))
        TRAJ.append(basemanip.MoveToHandPosition(matrices=T, outputtrajobj=True))
        waitrobot(robot)
    return planningutils.MergeTrajectories(TRAJ)      


def move_dijkstra(turbine, blade, organized_rays_list, interpolation):
    robot = turbine.robot
    joint_path_list = []
    time = interpolation/turbine.config.coating.coating_speed
    limits = robot.GetDOFVelocityLimits()*time
    deep = False
    
    for organized_rays in organized_rays_list:
        linear_interpolation = generate_linear_interpolation_rays(organized_rays, blade, interpolation)
        for deep in [False,True]:          
            try:
                joints = planning.joint_planning(turbine, linear_interpolation, deep)
            except IndexError:
                continue
            if len(joint_path_list)!=0:
                joints.insert(0,[joint_path_list[-1][-1]])
                joint_path, path, min_cost, adj, cost = planning.make_dijkstra(joints, limits, True)

                print min_cost
                if min_cost != inf:
                    joint_path_list.append(joint_path[1:])
                    break
##                joint_path_list.append(joint_path[1:]) # TIRA ISSO
                
            else:
                joint_path, path, min_cost, adj, cost = planning.make_dijkstra(joints, limits, True)
                print min_cost
                if min_cost != inf:
                    joint_path_list.append(joint_path)
                    break
##                joint_path_list.append(joint_path) # TIRA ISSO
        else:
            return [[]]
##            return joint_path_list
        
    return joint_path_list
                

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
