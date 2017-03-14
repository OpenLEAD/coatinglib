import cPickle
from turbine import Turbine
from turbine_config import TurbineConfig, ConfigFileError
import mathtools
from numpy import array, tan, arange
from math import pi as pi
from os.path import join, realpath
from os import environ
import itertools

def check_line(line, grids_num, grid_bases, line_grid, line_grid_dist):
    grid_dist = dict()
    for grid in grids_num:
        bases = list(grid_bases[grid])
        point_near, distance, distance_str = mathtools.distance_line_bases(
            line, bases, min_threshold, max_threshold)
        if distance!=None:
            line_grid[line] = line_grid.get(line,set()) | set([grid])
            grid_dist[grid] = [point_near, distance, distance_str]
    line_grid_dist[line] = grid_dist
    return line_grid, line_grid_dist

def compute_all_lines(a, b, grids_num, grid_bases):
    line_grid = dict()
    line_grid_dist = dict()
    for bi in b:
        for ai in a:
            line = (ai,bi)
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
                    print 'k=', k
                    print 'best_sol=', best_sol
                else:
                    best_sol.append(line_comb)
        if stop:
            break
    return best_sol          

if __name__ == '__main__':

    dir_test = join(realpath('.'),'test')
    environ['OPENRAVE_DATA'] = str(dir_test)
    cfg = TurbineConfig.load('turbine_unittest.cfg','test')
    turb = Turbine(cfg)

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
    a = tan(array(range(-80,80,1))*pi/180)
    b = arange(0,cfg.environment.y_min,-0.1)
    line_grid, line_grid_dist = compute_all_lines(a, b, grid_nums, grid_bases)
    best_sol = compute_minimal_lines(line_grid, grid_nums)
    
