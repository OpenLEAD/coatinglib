from numpy import array, linalg, dot, zeros
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
    organized_rays = [x for x in organized_rays if x != []]
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

def generate_linear_interpolation_joints(joint_solutions):
    new_joints = []
    new_joints.append(joint_solutions[0])
    for i in range(0,len(joint_solutions)-1):
        joints = mathtools.linear_interpolation_joint(joint_solutions[i], joint_solutions[i+1])
        new_joints.extend(joints[1:])
    return new_joints

def generate_linear_interpolation_rays(organized_rays):
    new_rays = []
    new_rays.append(organized_rays[0])
    for i in range(0,len(organized_rays)-1):
        points = mathtools.linear_interpolation_ray(organized_rays[i], organized_rays[i+1])
        new_rays.extend(points[1:])
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

def movetohandposition(robot, joint_solutions):
    manip = robot.GetActiveManipulator()
    basemanip = interfaces.BaseManipulation(robot,plannername='birrt')
    #TRAJ = []
    T = []
    for i in range(0,len(joint_solutions)-1):
        robot.SetDOFValues(joint_solutions[i])
        Ai = manip.GetTransform()
        robot.SetDOFValues(joint_solutions[i+1])
        Aj = manip.GetTransform()
        robot.SetDOFValues(joint_solutions[i])
        #Bi = logm(dot(linalg.inv(Ai),Aj))/3
        #Bj = logm(dot(linalg.inv(Aj),Ak))/3
        Bi = zeros((4,4))
        T.append(mathtools.homogenous_matrix_cubic_interpolation(Ai,Aj,Bi,Bi,10))
        #TRAJ.append(basemanip.MoveToHandPosition(
        #    matrices=mathtools.homogenous_matrix_cubic_interpolation(Ai,Aj,Bi,Bi,10), outputtrajobj=True))
        #waitrobot(robot)
    #traj = planningutils.MergeTrajectories(TRAJ)
    return T#traj
    

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
