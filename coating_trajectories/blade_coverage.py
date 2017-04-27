from numpy import array
import db
import planning

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

    DB.T = DB.turb.blades[0].GetTransform()
    parallels, borders = (DB.load_grid_to_trajectories())[grid]
    rays = DB.compute_rays_from_parallels(parallels, borders)
    blade = DB.load_blade()

    organized_rays = []
    for i in range(0,len(rays)):
        if i%2==0:
            organized_rays += rays[i]
        else:
            organized_rays += reversed(rays[i])
    organized_rays = [x for x in organized_rays if x != []]

    DB.turb.robot.GetLink('Flame').Enable(False)
    joint_solutions = planning.compute_robot_joints_opt(DB.turb, organized_rays, 0,
                                                    blade.trajectory_iter_surface)
    score = len(joint_solutions)*1.0/len(organized_rays)
    return score, joint_solutions

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
