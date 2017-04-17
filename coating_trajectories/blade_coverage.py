import cPickle
from turbine import Turbine
from turbine_config import TurbineConfig, ConfigFileError
import mathtools
from numpy import array, tan, arange, zeros, linalg, sign
from math import pi as pi
from os.path import join, realpath
from os import environ, makedirs
from visualizer import Visualizer
import db
import blade_modeling
import errno
import rail_place
from math import atan2
from itertools import combinations

def set_union_bases(line_grid, lines):
    line_union = set()
    for line in lines:
        line_union = line_union.union(line_grid[line])
    return line_union

def compute_minimal_lines(line_grid, grid_nums):
    lines = line_grid.keys()
    set_grid_nums = set(grid_nums)
    k = len(set_grid_nums)
    comb_sol = []
    best_sol = []
    stop=0
    for i in range(1,len(lines)):
        for line_comb in combinations(lines,i):
            line_union = set_union_bases(line_grid, line_comb)
            difference = set_grid_nums.difference(line_union)
            n = len(difference)
            if n==0:
                stop=1
            if n<=k:
                if n<k:
                    k = n
                    best_sol = []
                    best_sol.append(line_comb)
                else:
                    best_sol.append(line_comb)
        if stop:
            break
    return best_sol          

def remove_nonstr_lines(line_grid, line_grid_dist, threshold_str):
    for key in line_grid.keys():
        for grid in line_grid_dist[key].keys():
            [point_near, distance, distance_str] = line_grid_dist[key][grid]
            if distance_str < threshold_str:
                grids = line_grid[key]
                grids.remove(grid)
                line_grid[key] = grids
                line_grid_dist[key].pop(grid,None)
    return line_grid, line_grid_dist
                
def remove_grid_bases(line_grid, grid_nums):
    for line in line_grid.keys():
        for grid in line_grid[line]:
            grid_nums.remove(grid)
    return grid_nums

def plot_line(vis, lines, key):
    p0 = (-2,0,cfg.environment.z_floor_level)
    p1 = (2,0,cfg.environment.z_floor_level)
    points = array((p0,p1))
    vis.drawline(points, key, (0,0,0), 5)

    for line in lines:
        # Ax+By+C=0
        x1 = line[0][0]; y1 = line[0][1]
        x2 = line[1][0]; y2 = line[1][1]
        A = y1-y2; B = x2-x1
        C = x1*y2-x2*y1
        if A!=0:
            y1=cfg.environment.y_min; y2=cfg.environment.y_max
            p0 = ((-C-y1*B)/A,y1,cfg.environment.z_floor_level)
            p1 = ((-C-y2*B)/A,y2,cfg.environment.z_floor_level)
        else:
            p0 = (cfg.environment.x_min,y1,cfg.environment.z_floor_level)
            p1 = (cfg.environment.x_max,y1,cfg.environment.z_floor_level)
        points = array((p0,p1))
        vis.drawline(points, key, (0,0,1), 5)
    return

def plot_grids(vis, lines, line_grid):
    db_grid_to_trajectories =DB.load_db_grid_to_trajectories()
    blade = load_blade(blade_folder)
    grid_list = []
    line0 = (float('Inf'), float('Inf'))
    grids = list(line_grid[line0])
    grid_list.extend(grids)
    for line in lines:
        grid_list.extend(line_grid[line])
    for grid in grid_list:
        trajectories, borders = db_grid_to_trajectories[grid]
        rays = DB.compute_rays_from_parallels(blade, trajectories, borders)
        if grid_list.count(grid)>2:
            vis.plot_lists(rays,'rays',(0,0,1))
        if grid_list.count(grid)==1:
            vis.plot_lists(rays,'rays',(0,0,0))
        if grid_list.count(grid)==2:
            vis.plot_lists(rays,'rays',(1,0,0))
        grid_list = filter(lambda a: a != grid, grid_list)
    return 
        
def plot_bases(vis, sol, line_grid, line_grid_dist):
    line0 = (float('Inf'), float('Inf'))
    grids = line_grid[line0]
    for grid in grids:
        point_near, distance, distance_str = line_grid_dist[line0][grid]
        p = (point_near[0],point_near[1],cfg.environment.z_floor_level)
        vis.plot(p,'point',(1,0,0),10)
        
    for line in sol:
        grids = line_grid[line]
        for grid in grids:
            point_near, distance, distance_str = line_grid_dist[line][grid]
            p = (point_near[0],point_near[1],cfg.environment.z_floor_level)
            vis.plot(p,'point',(1,0,0),10)
    return 
            
def sort_best_sol(best_sol, line_grid, line_grid_dist):
    min_str = []
    mean_str = []
    sum_str = []
    num_grids = []
    
    for sol in best_sol:
        ngrids = 0
        distance_str_total = 0
        distance_str_min = 1000
        for line in sol:
            grids = line_grid[line]
            ngrids+=len(grids)
            for grid in grids:
                point_near, distance, distance_str = line_grid_dist[line][grid]
                distance_str_total+=distance_str
                distance_str_min = min(distance_str_min,distance_str)
        sum_str.append(distance_str_total)
        mean_str.append(distance_str_total*1.0/ngrids)
        num_grids.append(ngrids)
        min_str.append(distance_str_min)

    values = zeros((len(best_sol),4))
    values[:,0] = sum_str; values[:,1] = mean_str; values[:,2] = num_grids; values[:,3] = min_str
    best_sol_ordered = [x for (y,x) in sorted(zip(-values[:,0],best_sol))]
    return values, best_sol_ordered

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
    remove = [59]
    for i in remove:
        grid_nums.remove(i)
    return grid_nums

def lip_grids():
    return [0,1,2]

def jusante(x_range, grid_path):
    grid_path = 'jusante'
    grid_nums = jusante_grids()
    x_range = [0.7,2]
    return grid_nums, x_range, grid_path

def montante(x_range, grid_path):
    grid_path = 'montante'
    grid_nums = montante_grids()
    x_range = [-2,-0.5]
    return grid_nums, x_range, grid_path

def lip(x_range, grid_path):
    grid_path = 'lip'
    grid_nums = lip_grids()
    x_range = [0.7,2]
    return grid_nums, x_range, grid_path

def save_best_sol(best_sol_ordered, line_grid, line_grid_dist):
    path = join(db_directories,'rails',grid_path)
    try:
        makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise    
    DB.save_db_pickle(best_sol_ordered,join(path, 'best_sol.pkl'))
    DB.save_db_pickle(line_grid,join(path, 'line_grid.pkl'))
    DB.save_db_pickle(line_grid_dist,join(path, 'line_grid_dist.pkl'))            
    return

def robot_base_position(db_directories, grid_path, grid):
    path = join(db_directories,'rails',grid_path)
    DB = db.DB(db_directories)
    best_sol_ordered = DB.load_db_pickle(join(path, 'best_sol.pkl'))
    line_grid = DB.load_db_pickle(join(path, 'line_grid.pkl'))
    line_grid_dist = DB.load_db_pickle(join(path, 'line_grid_dist.pkl'))

    line0 = (float('Inf'), float('Inf'))
    grids = line_grid[line0]
    if grid in grids:
        point_near, distance, distance_str = line_grid_dist[line0][grid]
        return (point_near[0], 0, 0)
        
    for line in best_sol_ordered[0]:
        grids = line_grid[line]
        if grid in grids:
            # Ax+By+C=0
            x1 = line[0][0]; y1 = line[0][1]
            x2 = line[1][0]; y2 = line[1][1]
            point_near, distance, distance_str = line_grid_dist[line][grid]
            p = mathtools.closest_point_line_3d(array(line[0]), array(line[1]), point_near)
            return (x1, linalg.norm(p-line[0]), atan2(x1-x2,y2-y1)) 

def select_db(grids_to_coat):
    db_grids = DB.get_dbs_grids()
    for ncombination in range(1,len(db_grids)+1):
        feasible_combination = []
        for dbs in combinations(db_grids.keys(),ncombination):
            coatable = set([])
            for dbi in dbs:
                coatable = coatable|db_grids[dbi]
            for grid in grids_to_coat:
                if grid not in coatable:
                    break
            else: feasible_combination.append(dbs)
        if len(feasible_combination)>0:
            return feasible_combination
    return 
                
 
if __name__ == '__main__':

    dir_test = join(realpath('.'),'test')
    environ['OPENRAVE_DATA'] = str(dir_test)
    cfg = TurbineConfig.load('turbine_unittest.cfg','test')
    turb = Turbine(cfg)
    key='l'

    # Visualizer
    vis = Visualizer(turb.env)
    
    # DB inputs
    path = 'FACE'
    DB = db.DB(path)
    grids_to_coat = []

    with open(db_directories+'/grid_bases.pkl', 'rb') as f:
            grid_bases = cPickle.load(f)

    # Side inputs
    x_range = [0,0]
    angle_range = [-90,90] 
    #grid_nums, x_range, grid_path = jusante(x_range, grid_path)
    #grid_nums, x_range, grid_path = montante(x_range, grid_path)
    #grid_nums, x_range, grid_path = lip(x_range, grid_path)
    #lines = compute_lines(x_range, angle_range)

    # Line Parameters
    min_threshold, max_threshold = 0.1, 0.2
    
    # Primary rail computation
    threshold_str = 5
    #line_grid, line_grid_dist = primary_rail_grids(grid_nums, grid_bases)
    #line_grid, line_grid_dist = remove_nonstr_lines(line_grid, line_grid_dist, threshold_str)
    #grid_nums = remove_grid_bases(line_grid, grid_nums)
    
    # Secondary rail computation
    #line_grid, line_grid_dist = compute_all_lines(lines, grid_nums, grid_bases, line_grid, line_grid_dist)
    #best_sol = compute_minimal_lines(line_grid, grid_nums)
    #values, best_sol_ordered = sort_best_sol(best_sol, line_grid, line_grid_dist)
    #save_best_sol(best_sol_ordered, line_grid, line_grid_dist)
    #psa = robot_base_position(db_directories, grid_path, grid_num)
    #rp = rail_place.RailPlace(psa)
    #turb.place_rail(rp)

    
    
