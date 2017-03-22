from math import pi
import db
from os.path import join, realpath
import os
from turbine import Turbine
from turbine_config import TurbineConfig, ConfigFileError
from visualizer import Visualizer
from numpy import ones, array, dot, linspace
from openravepy import matrixFromAxisAngle
import rail_place
import mathtools
import blade_modeling
import planning
from copy import deepcopy
import csv
import cPickle

def generate_joints_db(turb, trajectory, blade):
    """
    Given the real blade angle:
    1) rotate the blades (update the environment);
    2) rotate the RBF model of the blade;
    3) rotate grid points;
    4) organize trajectories (removing empty points, adding borders,
    and making one full zigzagging list);
    5) compute optimization.

    Keyword arguments:
    blade_angle -- real angle of the blade
    rays -- points to be coated (list-n of lists-mx6) w.r.t. the 0 rotor angle.
    """
    

    turb.robot.GetLink('Flame').Enable(False)
    joint_solutions = planning.compute_robot_joints(turb, trajectory, 0, blade.trajectory_iter_surface)
    
    
    return joint_solutions


def generate_joints_db(trajectory, iter_surface, blade):
    """
    Method generate the coating trajectories. The trajectories are
    the intersection between two surfaces: the blade model, and the surface
    to be iterated, e.g. spheres. The algorithm
    follows the marching method, documentation available in:
    http://www.mathematik.tu-darmstadt.de/~ehartmann/cdgen0104.pdf
    page 94 intersection surface - surface.

    This is a data-intensive computing and might freeze your computer.

    Keyword arguments:
    iter_surface -- surface to be iterated, as mathtools.sphere.
    step -- it must be small, e.g. 1e-3. Otherwise the method will fail.
    """

    step = self.trajectory_step

    trajectories = deepcopy(self.trajectories)
    
    if not issubclass(iter_surface.__class__, mathtools.IterSurface):
        raise TypeError("Object is not a valid surface.")

    if self.model_iter_surface is not None:
        if iter_surface.name() != self.model_iter_surface.name():
            raise TypeError("Object iter_surface must have same type of model_iter_surface.")

    if not self.models:
        raise IndexError("Object models is empty. Load or create a model before generate trajectories.")

    self._blade.SetTransform(eye(4))
    self._blade.SetTransform(dot(self._blade.GetTransform(),
                                 matrixFromAxisAngle([0, -self.turbine.config.environment.blade_angle, 0]))
                             )

    for point in trajectory:
        model = blade.select_model(point_on_surfaces)
        
        
    point_on_surfaces = self.compute_initial_point(iter_surface, trajectories)
    self.trajectory_iter_surface = iter_surface

    while iter_surface.criteria():
        model = self.select_model(point_on_surfaces)
        trajectories.append(self.draw_parallel(point_on_surfaces, model, iter_surface, step))
        p0=trajectories[-1][-1]
        iter_surface.update()
        point_on_surfaces = mathtools.curvepoint(model, iter_surface, p0[0:3])
        self.trajectories = trajectories
    return trajectories



if __name__ == '__main__':

    #-----------------------------------------------------------------------
    # DB inputs
    db_directories = ['db']#, 'db_45', 'db_-45']
    db_angles = [0]#, pi/4,-pi/4]
    blade_folder = ['jiraublade_hd_filtered']#, 'jiraublade_hd_filtered_45',
                    #'jiraublade_hd_filtered_-45']
    #-----------------------------------------------------------------------
    
    dir_test = join(realpath('.'),'test')
    os.environ['OPENRAVE_DATA'] = str(dir_test)
    cfg = TurbineConfig.load('turbine_unittest.cfg','test')
    turb = Turbine(cfg)

    dict_angle_db = dict()
    dict_angle_blade = dict()
    for i in range(0,len(db_directories)):
        dict_angle_db[db_angles[i]]=db.DB(db_directories[i])
    for i in range(0,len(blade_folder)):
        dict_angle_blade[db_angles[i]]=load_blade(blade_folder[i])
    
    #vis = Visualizer(turb.env)

    # Grid input
    #grid_num = 26
    #sorted_bases, trajectories, borders = base_for_grid_coating(grid_num)
    #plot_base_for_grid_coating(sorted_bases, trajectories, borders)
    while True:
        if not tolerance_test(): break
    #non_coatable_one_base, non_coatable = non_coatable_grids()
    #non_coatable = plot_non_coatable_points_grid(grid_num)
