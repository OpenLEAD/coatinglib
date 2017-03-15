import cPickle
from turbine import Turbine
from turbine_config import TurbineConfig, ConfigFileError
import mathtools
from numpy import array, tan, arange
from math import pi as pi
from os.path import join, realpath
from os import environ
import itertools
from visualizer import Visualizer

def check_line(line, grids_num, grid_bases, line_grid, line_grid_dist):
    x1 = line[0]; x2 = line[1]
    grid_dist = dict()
    for grid in grids_num:
        bases = list(grid_bases[grid])
        point_near, distance, distance_str = mathtools.distance_line_bases(
            x1, x2, bases, min_threshold, max_threshold)
        if distance!=None:
            line_grid[line] = line_grid.get(line,set()) | set([grid])
            grid_dist[grid] = [point_near, distance, distance_str]
    line_grid_dist[line] = grid_dist
    return line_grid, line_grid_dist

def compute_all_lines(lines, grids_num, grid_bases, line_grid = dict(), line_grid_dist = dict()):
    for line in lines:
        line_grid, line_grid_dist = check_line(line, grids_num,
                                          grid_bases, line_grid, line_grid_dist)
    return line_grid, line_grid_dist

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
        for line_comb in itertools.combinations(lines,i):
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

def primary_rail_grids(grids_num, grid_bases, line_grid=dict(), line_grid_dist=dict()):
    grid_dist = dict()
    x1 = (0,0)
    x2 = (1,0)
    line = (float('Inf'), float('Inf'))
    
    for grid in grids_num:
        bases = list(grid_bases[grid])
        point_near, distance, distance_str = mathtools.distance_line_bases(
            x1, x2, bases, min_threshold, max_threshold)
        if distance!=None:
            line_grid[line] = line_grid.get(line,set()) | set([grid])
            grid_dist[grid] = [point_near, distance, distance_str]
    line_grid_dist[line] = grid_dist
    return line_grid, line_grid_dist

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
    p0 = (0,0,cfg.environment.z_floor_level)
    p1 = (2,0,cfg.environment.z_floor_level)
    points = array((p0,p1))
    vis.drawline(points, key, (0,0,0), 5)
    for line in lines:
        p0 = (line[0][0],line[0][1],cfg.environment.z_floor_level)
        p1 = (line[1][0],line[1][1],cfg.environment.z_floor_level)
        points = array((p0,p1))
        vis.drawline(points, key, (0,0,1), 5)
    return
    
def compute_lines(x_range, angle_range):
    lines=[]
    for x in arange(x_range[0],x_range[1],0.1):
        for angle in arange(angle_range[0],-8,8):
            y = tan(angle*pi/180)
            if y < cfg.environment.y_max:
                lines.append(((x,0),(0,-x*y)))
        for angle in arange(8,angle_range[1],8):
            y = tan(angle*pi/180)
            if y > cfg.environment.y_min:
                lines.append(((x,0),(0,-x*y)))
        lines.append(((x,0),(x,1)))
    return lines

if __name__ == '__main__':

    dir_test = join(realpath('.'),'test')
    environ['OPENRAVE_DATA'] = str(dir_test)
    cfg = TurbineConfig.load('turbine_unittest.cfg','test')
    turb = Turbine(cfg)
    key='l'

    with open('grid_bases.pkl', 'rb') as f:
            grid_bases = cPickle.load(f)
    grid_nums = range(0,15)
    grid_nums.extend(range(17,20))
    grid_nums.extend(range(22,24))
    grid_nums.extend(range(67,69))
    grid_nums.extend(range(72,77))
    grid_nums.extend(range(78,80))
    grid_nums.remove(74)

    min_threshold, max_threshold = 0.1, 0.2
    threshold_str = 5
    x_range = [0.7,2]
    angle_range = [-80,80] 
    lines = compute_lines(x_range, angle_range)
    
    line_grid, line_grid_dist = primary_rail_grids(grid_nums, grid_bases)
    line_grid, line_grid_dist = remove_nonstr_lines(line_grid, line_grid_dist, threshold_str)
    grid_nums = remove_grid_bases(line_grid, grid_nums)
    
    line_grid, line_grid_dist = compute_all_lines(lines, grid_nums, grid_bases, line_grid, line_grid_dist)
    best_sol = compute_minimal_lines(line_grid, grid_nums)

    # Visualizer
    vis = Visualizer(turb.env)
    
