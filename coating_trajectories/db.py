#!/usr/bin/env python
from path_filters import filter_points
import robot_utils
from numpy import random, array, linspace, cross, transpose
from numpy import sign, dot, linalg, sum, zeros, round, delete
from os.path import splitext, join, exists, isfile, split
from os import makedirs, listdir
import copy
import rail_place
import cPickle
import mathtools
import errno
from lxml import etree as ET
from openravepy import matrixFromAxisAngle
import blade_modeling
from scipy.spatial import ConvexHull
from itertools import combinations
from math import atan2, pi

## @file
# @brief This contains functions to manage the database.
# @author Renan S. Freitas
# @bug No known bugs


class NoDBFound(Exception):
    def __init__(self, value):
        Exception.__init__(self, "No DB file named " + value + " found.")


def save_pickle(obj, name):
    """ Save a general file in pickle format.

    Args:
        obj: object to be saved.
        name: name of the file ('example/file.pkl').

    Examples:
        >>> save_pickle(my_var, 'my_var.pkl')
    """

    with open(name, 'wb') as f:
        cPickle.dump(obj, f, cPickle.HIGHEST_PROTOCOL)
    return


def load_pickle(filename):
    """ Load a general file in pickle format.

    Args:
        filename: path/name of the file.

    Examples:
        >>> load_pickle('my_folder/my_file.pkl')
    """

    with open(filename, 'rb') as f:
        return cPickle.load(f)


def get_sorted(filename):
    """ Given the base_to_num or point_to_num files, the method returns the sorted list, i.e., num_to_base or
    num_to_point.

    Args:
        filename: path/name of the file.

    Examples:
        >>> get_sorted('my_folder/base_to_num.pkl')
    """
    b2n = load_pickle(filename)
    return [b for (v, b) in sorted(zip(b2n.values(), b2n.keys()))]


class DB:
    """ DB for robot's base position and coated points.
    @page DB

    Args:
        path: is the folder to store all data bases;
        turbine: (@ref Turbine) turbine object;
        db_name: (optional) name of the db if it already exists.

    Examples:
        >>> DB = db.DB('FACE', turbine, 'db')
    """

    def __init__(self, path, turbine, db_name=None):

        self.path = path
        self.db_main_path = ''
        self.T = matrixFromAxisAngle([0, 0, 0])
        self.info = None
        self.turbine = turbine
        self.verticalized = False

        if not exists(self.path):
            makedirs(self.path)

        try:
            self.info = ET.parse(open(join(self.path, 'info.xml')))
        except IOError as error:
            if error.errno == errno.ENOENT:
                raise NoDBFound('info.xml')
        self.verticalization_state()

        if db_name != None and db_name != '':
            dbs = self.info.findall('db')
            for db in dbs:
                if db.find('name').text == db_name:
                    self.db_main_path = join(self.path, db.find('path').text)
                    self.T = self._extract_T(db)
                    break
            else:
                raise NoDBFound(db_name)
        else:
            self.db_main_path = self.path

    def verticalization_state(self):
        """ Check if db was verticalized. DB.verticalized == True if verticalized base.
        The Lip and borders are, for example, verticalized DBs.

        Examples:
            >>> DB.verticalized()
        """

        vertical = self.info.find('verticalized')
        if vertical is None:
            self.verticalized = False
            return
        if vertical.text == '1':
            self.verticalized = True
            return

    def _extract_T(self, info_db):
        """ Return the homogeneous matrix of the blade when the base was computed (rotation and position).

         Args:
             info_db: xml file of the DB.

        Examples:
            >>> info = ET.parse(open(join(self.path, 'info.xml')))
            >>> m = DB._extract_T(info)
        """

        T = info_db.find('transform').text
        T = T.replace('\n', ' ')
        T = T.replace('\t', ' ')
        T = T.split(' ')
        T[:] = [x for x in T if x != '']
        m = []
        mi = []
        for t in T:
            try:
                mi.append(float(t))
            except:
                None
            if len(mi) == 4:
                m.append(mi)
                mi = []
        m = array(m)
        if m.shape == (4, 4):
            return m
        else:
            raise SyntaxError("Invalid transform shape in info.xml")

    def load_db(self):
        """ Load main database num:num. The first num represents the points and
        the second represents the bases.

        Examples:
            >>> DB.load_db()
        """

        name = join(self.db_main_path, 'db.pkl')
        return load_pickle(name)

    def load_points_to_num(self):
        """ Load points_to_num database points:num. Points are tuples (x,y,z) and
        the num maps the points (counter to reduce database complexity). It will rotate and
        translate the points by the T matrix (database homogeneous transform).

        Examples:
            >>> DB.load_points_to_num()
        """

        if self.verticalized:
            path = join(self.path, 'points_to_num_verticalized.pkl')
        else:
            path = join(self.path, 'points_to_num.pkl')
        ptn = load_pickle(path)
        real_ptn = dict()
        for key, value in ptn.iteritems():
            real_ptn[tuple(round(dot(self.T, [key[0], key[1], key[2], 1])[0:3], 9))] = value
        return real_ptn

    def load_points_to_num_standard(self):
        """ Load points_to_num database points:num. Points are tuples (x,y,z) and
        the num maps the points (counter to reduce database complexity). It will rotate and
        translate the points by the T matrix (database homogeneous transform).

        Examples:
            >>> DB.load_points_to_num()
        """

        path = join(self.path, 'points_to_num.pkl')
        ptn = load_pickle(path)
        real_ptn = dict()
        for key, value in ptn.iteritems():
            real_ptn[tuple(round(dot(self.T, [key[0], key[1], key[2], 1])[0:3], 9))] = value
        return real_ptn

    def load_grid_to_mp(self):
        """ Load grid_to_mp database num:[(m1,m2),(p1,p2)].
        Grids are numbers and the mp are
        meridian/parallels index.

        Examples:
            >>> DB.load_grid_to_mp()
        """

        path = join(self.path, 'grid_to_mp.pkl')
        return load_pickle(path)

    def load_grid_to_trajectories(self):
        """ Load grid_to_trajectories database num:[trajectories, borders]. Grids are numbers and the bases are lists
        of trajectories (num, db format) and borders (tuples Nx3). It will rotate and translate the points by the
        T matrix (database homogeneous transform). It will check if the trajectories are verticalized and it will load
        the correct file.

        Examples:
            >>> DB.load_grid_to_trajectories()
        """

        if self.verticalized == True:
            path = join(self.path, 'grid_to_trajectories_verticalized.pkl')
        else:
            path = join(self.path, 'grid_to_trajectories.pkl')
        grid_to_trajectories = load_pickle(path)
        new_grid_to_trajectories = dict()
        for grid, value in grid_to_trajectories.iteritems():
            trajectories, borders = value
            new_grid_to_trajectories[grid] = [trajectories, mathtools.rotate_points(borders, self.T)]
        return new_grid_to_trajectories

    def load_grid_to_trajectories_standard(self):
        path = join(self.path, 'grid_to_trajectories.pkl')
        grid_to_trajectories = load_pickle(path)
        new_grid_to_trajectories = dict()
        for grid, value in grid_to_trajectories.iteritems():
            trajectories, borders = value
            new_grid_to_trajectories[grid] = [trajectories, mathtools.rotate_points(borders, self.T)]
        return new_grid_to_trajectories

    def load_bases_to_num(self):
        """ Load bases_to_num database bases:num. Bases are tuples railplace and
        the num maps the bases (counter to reduce database complexity).

        Examples:
            >>> DB.load_bases_to_num()
        """

        path = join(self.path, 'bases_to_num.pkl')

        return load_pickle(path)

    def load_visited_bases(self):
        """ Load visited_bases database bases:bool. Bases are tuples PSAlpha and
        the bool show if it was already computed.

        Examples:
            >>> DB.load_visited_bases()
        """

        path = join(self.db_main_path, 'visited_bases.pkl')
        return load_pickle(path)

    def load_grid_meridian(self):
        """ Load grid_meridians. The meridians used to make the grid.

        Examples:
            >>> DB.load_grid_meridian()
        """

        path = join(self.path, 'meridians.pkl')
        return mathtools.rotate_trajectories(load_pickle(path), self.T)

    def load_grid_parallel(self):
        """ Load grid_parallel. The parallels used to make the grid.

        Examples:
            >>> DB.load_grid_parallel()
        """

        path = join(self.path, 'parallels.pkl')
        return mathtools.rotate_trajectories(load_pickle(path), self.T)

    def load_blade(self):
        """ Function to load blade model.

        Examples:
            >>> DB.load_blade()
        """

        folder = self.info.find('blade_model_filtered')
        if folder == None:
            folder = self.info.find('blade_model')
        if folder == None:
            raise NameError("No blade model found in info.xml")

        folder = folder.text
        xml_trajectories_path = join(self.path, folder, "trajectory/trajectory.xml")
        blade = blade_modeling.BladeModeling(self.turbine, self.turbine.blades[3])
        blade.load_trajectory(xml_trajectories_path)
        blade.trajectories = mathtools.rotate_trajectories(blade.trajectories, self.T)
        blade.rotate_models(self.T)
        return blade

    def load_blade_full(self):
        """ Function to load blade model.

        Examples:
            >>> DB.load_blade_full()
        """

        folder = self.info.find('blade_model')
        if folder == None:
            raise NameError("No blade model found in info.xml")

        folder = folder.text
        xml_trajectories_path = join(self.path, folder, "trajectory/trajectory.xml")
        blade = blade_modeling.BladeModeling(self.turbine, self.turbine.blades[0])
        blade.load_trajectory(xml_trajectories_path)
        blade.trajectories = mathtools.rotate_trajectories(blade.trajectories, self.T)
        blade.rotate_models(self.T)
        return blade

    def create_points_to_num(self):
        """ This method get the samples from the blade and create the point_to_num database: a dictionary
        that relates samples (keys) and numbers (values, counters to minimize the database complexity).
        It will return the point_to_num dictionary (it will not save it).

        Examples:
            >>> DB.create_points_to_num()
        """

        blade = self.load_blade()
        points_to_num = dict()

        counter = 0
        for trajectory in blade.trajectories:
            for point in trajectory:
                points_to_num[tuple(round(point[0:3], 9))] = counter
                counter += 1
        return points_to_num

    def create_db(self):
        """ Create the database.
        If point_to_num does not exist, it will create it and save it.
        It will make a single db looping on all databases (angles)

        Examples:
            >>> DB.create_db()
        """

        try:
            points_to_num = self.load_points_to_num()
        except IOError:
            points_to_num = self.create_points_to_num()
            save_pickle(points_to_num, join(self.path, 'points_to_num.pkl'))

        try:
            main_db = self.load_db()
        except IOError:
            main_db = dict()
            if self.path == self.db_main_path:
                dbs = self.info.findall('db')
                for dbi in dbs:
                    T = self._extract_T(dbi)
                    db_path = join(self.path, dbi.find('path').text)
                    main_dbi = load_pickle(join(db_path, 'db.pkl'))
                    for point, base in main_dbi.iteritems():
                        b = main_db.get(point, [])
                        b.append([base, T])
                        main_db[point] = b
            else:
                for k,v in points_to_num.iteritems():
                    main_db[v] = set([])
            save_pickle(main_db, join(self.db_main_path, 'db.pkl'))
        del main_db

        try:
            _ = self.load_bases_to_num()
        except IOError as error:
            if error.errno == errno.ENOENT:
                raise NoDBFound('bases_to_num.pkl')
            else:
                raise
        return

    def create_grid_verticalization(self, n=60, distance=3e-3):
        """ This method iterates on grids. For each grid, it creates interpolations (with legMLS) for all parallels,
        making more points for verticalization.

        Args:
            grid_to_trajectories: (dict) key: grids, values: [trajectories, borders]
            points_to_num: (dict) key: (tuple) points, values: (int) num
            n: number of points per trajectory.
            distance: distance between points.
        """

        blade = self.load_blade()
        counter = 0
        grid_to_trajectories = self.load_grid_to_trajectories_standard()
        ptn = self.load_points_to_num_standard()
        points_to_num_verticalized = dict()
        for grid, value in grid_to_trajectories.iteritems():
            traj, bord = value
            rays_list = self.compute_rays_from_parallels(traj, bord, ptn)
            parallels = []
            for r, rays in enumerate(rays_list):
                if len(rays) <= 3:
                    continue
                lin = linspace(0, len(rays) - 1, n)
                N = len(lin)
                new_rays = zeros((N, 6))
                new_rays[:, 0], _, _ = mathtools.legMLS(
                    array(rays)[:, 0], array(range(len(rays))), lin, 3, scale=1.5, dwf=None, ddwf=None)
                new_rays[:, 1], _, _ = mathtools.legMLS(
                    array(rays)[:, 1], array(range(len(rays))), lin, 3, scale=1.5, dwf=None, ddwf=None)
                new_rays[:, 2], _, _ = mathtools.legMLS(
                    array(rays)[:, 2], array(range(len(rays))), lin, 3, scale=1.5, dwf=None, ddwf=None)
                remove = []
                for i in range(len(new_rays)):
                    try:
                        new_rays[i] = blade.compute_ray_from_point(new_rays[i][0:3])
                    except ValueError:
                        remove.append(i)
                        continue
                new_rays = delete(new_rays, remove, 0)
                parallels.append(new_rays)
            parallels = mathtools.equally_spacer(parallels, distance=distance)
            parallels = mathtools.trajectory_verticalization(parallels)

            R = transpose(self.T[0:3, 0:3])
            borders = []
            for i in range(len(parallels)):
                borders.append([dot(R,parallels[i][0][0:3]), dot(R,parallels[i][-1][0:3])])
                parallels[i] = parallels[i][1:-1]

            parallels_num = []
            for parallel in parallels:
                parallel_num = []
                for point in parallel:
                    point[0:3] = dot(R,point[0:3])
                    points_to_num_verticalized[tuple(round(point[0:3], 9))] = counter
                    parallel_num.append(counter)
                    counter += 1
                parallels_num.append(parallel_num)

            grid_to_trajectories[grid] = [parallels_num, borders]

        return grid_to_trajectories, points_to_num_verticalized

    def get_sorted_bases(self):
        """ Return the sorted bases -- tuples PSAlpha.

        Examples:
            >>> ntb = DB.get_sorted_bases()
        """

        btn = self.load_bases_to_num()
        return [b for (v, b) in sorted(zip(btn.values(), btn.keys()))]

    def get_sorted_points(self, points_to_num = None):
        """ Return the sorted points -- tuples (x,y,z).

        Examples:
            >>> ntp = DB.get_sorted_points()
        """

        if points_to_num is None:
            points_to_num = self.load_points_to_num()
        return [k for (v, k) in sorted(zip(points_to_num.values(), points_to_num.keys()))]

    @staticmethod
    def invert_dict(dict_a):
        dict_b = dict()
        for k, v in dict_a.iteritems():
            dict_b[v] = k
        return dict_b


    def get_bases(self, db):
        """ Return all feasible bases in given database

        Examples:
            >>> bases = DB.get_bases(db)
        """

        bases = set()
        for val in db.values():
            bases = val | bases
        return bases

    def base_grid_coating(self, grid_num, grid_to_trajectories = None):
        """ Method returns the bases in score order that can coat a specific grid.

        Args:
            grid_num: grid to be coated

        Examples:
            >>> bases, scores = DB.base_grid_coating(grid)
        """

        if grid_to_trajectories is None:
            grid_to_trajectories = self.load_grid_to_trajectories()

        trajectories, borders = grid_to_trajectories[grid_num]
        bases, scores = self.base_trajectory_coating(trajectories)
        return bases, scores

    def base_trajectory_coating(self, trajectories):
        """ Method returns the bases in score order (0 to 1) that can coat the given
        trajectories. This order is more important to know than absolut coating
        because some points may never be coated (as points on the top of
        the blade)

        Args:
            trajectories: num in points_to_num (db format) to be coated

        Examples:
            >>> best_bases, scores = DB.base_trajectory_coating(trajectories)
        """

        db = self.load_db()
        bases_tuple = self.get_sorted_bases()
        score = zeros(len(bases_tuple))
        N = 0

        for trajectory in trajectories:
            for point in trajectory:
                try:
                    score[list(db[point])] -= 1
                    N += 1
                except KeyError:
                    pass
        best_bases = [x for (y, x) in sorted(zip(score, bases_tuple))]
        return best_bases, -array(sorted(score)) * 1.0 / N

    def create_db_grid_bases(self, threshold):
        """ Creates a dictionary grid_num:bases which relates grid and feasible bases that can coat it given a
        threshold for the score.

        Args:
             threshold: (int) threshold for the score

        Examples:
            >>> DB.create_db_grid_bases(0.9)
        """

        gtt = self.load_grid_to_trajectories()
        grid_bases = dict()
        for grid_num in gtt.keys():
            bases, scores = self.base_grid_coating(grid_num, gtt)
            for i, base in enumerate(bases):
                if scores[i] >= threshold:
                    rp = rail_place.RailPlace(base)
                    xyz = rp.getXYZ()
                    value = {(xyz[0], xyz[1])}
                    grid_bases[grid_num] = grid_bases.get(grid_num, set()) | value
                else:
                    break
        save_pickle(grid_bases, join(self.db_main_path, 'grid_bases.pkl'))
        return

    def get_dbs_grids(self):
        """ Load grid_bases.pkl and info.xml files to return all coatable grids. The method returns a dictionary
        {db_path : coatable_grids}.

        Examples:
            >>> db_grids = DB.get_dbs_grids()
        """

        dbs = self.info.findall('db')
        db_grids = dict()
        for dbi in dbs:
            db_path = join(self.path, dbi.find('path').text)
            grid_bases = load_pickle(join(db_path, 'grid_bases.pkl'))
            grids = set([])
            for grid, base in grid_bases.iteritems():
                if len(base) > 0:
                    grids = grids | {grid}
            db_grids[db_path] = grids
        return db_grids

    def merge_db(self, db_file, main_db):
        """ Merge new db_file to the main database.

        Args:
            db_file: file to be merged;
            main_db: main database file;

        Examples:
            >>> main_db = DB.merge_db(db_file, main_db)
        """

        for key, value in db_file.iteritems():
            if type(value) == int:
                value = set([value])
            main_db[key] = main_db.get(key, set()) | value
        return main_db

    def merge_seg(self, seg, main_db):
        """ Merge new segments to the main database. Segments are parts of the parallels.

        Args:
            seg: segments to be merged
            main_db: main database file

        Examples:
           >>> main_db = DB.merge_seg(seg, main_db)
        """
        if isinstance(seg,dict):
            base = seg.keys()[0]
            for segments in seg[base]:
                for segment in segments:
                    for point in segment:
                        main_db[point] = main_db.get(point, set()) | set([base])
            return main_db
        else:
            return main_db

    def create_db_from_segments(self, path, merge=0):
        """ Method to create db (num to num) with segments in given path. Segments are parts of the parallels.
        It returns the resulted main_db (the method does not save it)

        Args:
            path: where the segments are.
            merge: (default=0) create new empty db if merge == 0. Else, it merges.

        Examples:
            >>> db = DB.create_db_from_segments(path)
        """

        if merge == 0:
            db = dict()
            ptn = self.load_points_to_num()
            for key, value in ptn.iteritems():
                db[value] = set()
        else:
            db = self.load_db()

        onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
        for afile in onlyfiles:
            filename, file_extension = splitext(afile)
            if file_extension == '.pkl':
                try:
                    seg = load_pickle(join(path, afile))
                except EOFError:
                    continue
                db = self.merge_seg(seg, db)
        return db

    def compile_db(self, path, merge=0):
        """ Method to compile dbs (num to num) generated by generate_db.
        It returns the resulted main_db.

        Args:
            path: where the dbs are.
            merge: (default=0) create new empty db if merge == 0. Else, the method will merge it with previous db.

        Examples:
            >>> db_main = DB.compile_dbs(path)
        """

        if merge == 0:
            db = dict()
            ptn = self.load_points_to_num()
            for key, value in ptn.iteritems():
                db[value] = set()
        else:
            db = self.load_db()

        onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
        for afile in onlyfiles:
            filename, file_extension = splitext(afile)
            if file_extension == '.pkl':
                try:
                    seg = load_pickle(join(path, afile))
                except EOFError:
                    continue
                db = self.merge_db(seg, db)
        return db

    def get_bases_trajectories(self, trajectories):
        """ Method returns the robot bases (tuples PSAlpha) that can coat all set of points in trajectories.
        Trajectories here are natural numbers (num in point_to_num).
        If a point in trajectories does not belong to the db, it is skipped.

        Args:
            trajectories: points_num (db format) to be coated

        Examples:
            >>> bases = DB.get_bases_trajectories(trajectories)
        """

        db = self.load_db()
        bases_tuple = self.get_sorted_bases()

        set_of_feasible_bases_num = self.get_bases(db)

        for trajectory in trajectories:
            for point in trajectory:
                try:
                    set_of_feasible_bases_num &= db[point]
                except KeyError:
                    pass
        return array(bases_tuple)[list(set_of_feasible_bases_num)].tolist()

    def make_grid(self, number_of_meridians, number_of_parallels, init_parallel):
        """ Make a grid in the blade with parallels and meridians.
        Meridians and parallels are evenly spaced.
        The parallel used to compute the meridians is in the middle.
        
        Args:
            number_of_meridians: number of meridians to be computed.
            number_of_parallels: number of parallels to be computed.
            init_parallel: initial parallel. This argument is important to remove the lip from the grid.

        Examples:
            >>> meridians, parallels = DB.make_grid(10, 10, 5)
        """

        blade = self.load_blade_full()
        parallel = blade.trajectories[int(len(blade.trajectories) / 2)]
        meridians = blade.draw_meridians(parallel, 1e-3, number_of_meridians)

        parallels = []
        list_index = linspace(init_parallel, len(blade.trajectories), number_of_parallels).astype(int)
        list_index[-1] -= 1
        for i in list_index:
            parallels.append(blade.trajectories[i])
        return meridians, parallels

    def create_db_grid(self):
        """ Autonomously create the db_grid after the make_grid process.
        The method will get adjacent meridians and parallels, building the grids.
        Some of the grids may not make sense, thus the user should analyse, add and remove grids.

        Examples:
            >>> DB.create_db_grid()
        """

        blade = self.load_blade_full()
        grid_to_mp = dict()
        grid_to_trajectories = dict()
        meridians = self.load_grid_meridian()
        parallels = self.load_grid_parallel()

        counter = 0
        for i in range(0, len(meridians)):
            for j in range(0, len(parallels) - 1):
                grid = [(i, (i + 1) % len(meridians)), (j, j + 1)]
                trajectories_in_grid, borders = self.get_points_in_grid(
                    blade, [meridians[grid[0][0]], meridians[grid[0][1]]],
                    [parallels[grid[1][0]], parallels[grid[1][1]]])

                grid_to_mp[counter] = grid
                grid_to_trajectories[counter] = [trajectories_in_grid, borders]
                counter += 1

        try:
            save_pickle(grid_to_mp, join(self.path, 'grid_to_mp.pkl'))
        except IOError:
            None

        try:
            save_pickle(grid_to_trajectories, join(self.path, 'grid_to_trajectories.pkl'))
        except IOError:
            None

        return

    def generate_db(self, base_num, points = None, do_side_filter = True):
        """ This method generates the database. Given the base (natural number that corresponds to a base
        tuple PSAlpha), it will test all points. It verifies only kinematic and collisions. It returns None if
        there is collision of the robot base and environment.

        It returns a dictionary ptb = point:base

        Warning: HARD COMPUTATION time

        Args:
            base_num: number in number_to_base
            do_side_filter: (optional, default = True) it will only check feasible segments on the side of the robot.
            It is a way of make things faster, but it will not work with borders and lip
        """

        ntb = self.get_sorted_bases()
        base = ntb[base_num]
        rp = rail_place.RailPlace(base)
        self.turbine.place_rail(rp)
        self.turbine.place_robot(rp)
        if self.turbine.check_rail_collision():
            return None
        if self.turbine.check_robotbase_collision():
            return None

        ptn = self.load_points_to_num()
        blade = self.load_blade()
        ptb = dict()

        if points is None:
            points = map(lambda x: blade.compute_ray_from_point(x), ptn.keys())

        points = filter_points(self.turbine, points, do_side_filter)

        for point in points:
            iksol = robot_utils.ik_angle_tolerance(self.turbine, point, deep = False)
            if len(iksol) > 0:
                ptb[ptn[tuple(round(point[0:3],9))]] = set([base_num])
        return ptb


    def get_points_in_grid(self, meridian, parallel, full=False):
        """ Get parallels that belong to a grid, between given meridians and given parallels.
        
        Args:
            meridian: tuple 1x2, First meridian must be on the right w.r.t. the second meridian.
            It means that points[1] (y) of points meridian1 > points[1] (y) of points meridian2.
            parallel: tuple 1x2. Not ordered.
        """

        if full:
            blade = self.load_blade_full()
        else:
            blade = self.load_blade()

        meridian1, meridian2 = meridian[0], meridian[1]
        parallel1, parallel2 = parallel[0], parallel[1]
        points_to_num = self.load_points_to_num()

        parallel_index_1 = 0
        parallel_index_2 = 0
        for i in range(0, len(blade.trajectories)):
            traj_list = [list(a) for a in blade.trajectories[i]]
            if list(parallel1[0]) in traj_list:
                parallel_index_1 = i
                break
        for i in range(0, len(blade.trajectories)):
            traj_list = [list(a) for a in blade.trajectories[i]]
            if list(parallel2[0]) in traj_list:
                parallel_index_2 = i
                break

        init = min(parallel_index_1, parallel_index_2)
        end = max(parallel_index_1, parallel_index_2) + 1
        trajectories_in_grid = []

        def get_point_value(point):
            try:
                return points_to_num[round(point, 9)]
            except KeyError:
                return None

        borders = []
        for i in range(init, end):
            trajectory_in_grid = []
            parallel = array(blade.trajectories[i])
            p1, sorted_parallel1 = self._closest_meridian_point(meridian1, parallel, blade)
            p2, sorted_parallel2 = self._closest_meridian_point(meridian2, parallel, blade)
            parallel1 = mathtools.direction_in_halfplane(parallel, p1[3:6])
            parallel2 = mathtools.direction_in_halfplane(parallel, p2[3:6])
            p1, sorted_parallel1 = self._closest_meridian_point(meridian1, parallel1, blade)
            p2, sorted_parallel2 = self._closest_meridian_point(meridian2, parallel2, blade)
            index_left = self._get_first_left_meridian_point_index(parallel, sorted_parallel1, p1)
            index_right = self._get_first_right_meridian_point_index(parallel, sorted_parallel2, p2)

            if index_left == None or index_right == None:
                pass
            elif abs(index_right - index_left) % (len(parallel) - 1) == 1:
                pass
            elif index_left <= index_right:
                trajectory_in_grid += filter(lambda x: x is not None,
                                             map(get_point_value,
                                                 map(tuple, array(parallel)[index_left:index_right + 1, 0:3].tolist())))
            else:
                trajectory_in_grid += filter(lambda x: x is not None,
                                             map(get_point_value,
                                                 map(tuple, array(parallel)[index_left:, 0:3].tolist()))) + \
                                      filter(lambda x: x is not None,
                                             map(get_point_value,
                                                 map(tuple, array(parallel)[:index_right + 1, 0:3].tolist())))
            trajectories_in_grid.append(trajectory_in_grid)
            borders.append([tuple(p1[0:3]), tuple(p2[0:3])])
        return trajectories_in_grid, borders

    def _closest_meridian_point(self, meridian, parallel, blade):
        min_dist = 100
        closest_meridian_point = []
        sorted_parallel = []
        dist_list = []
        for meridian_point in meridian:
            dist = sum((parallel[:, 0:3] - meridian_point[0:3]) *
                       (parallel[:, 0:3] - meridian_point[0:3]), 1)
            if min(dist) <= min_dist:
                closest_meridian_point = meridian_point
                min_dist = min(dist)
                dist_list = dist
        sorted_parallel = [x for (y, x) in sorted(zip(dist_list, parallel))]

        model = blade.select_model(closest_meridian_point)

        return blade.compute_ray_from_point(closest_meridian_point, model), sorted_parallel

    def _get_first_left_meridian_point_index(self, parallel, sorted_parallel, meridian_point):
        """
        This method assumes the use of spheres as iter_surface.
        """

        tan = cross(meridian_point[3:6],
                    meridian_point[0:3] / linalg.norm(meridian_point[0:3]))
        for point in sorted_parallel:
            if sign(dot(tan, meridian_point[0:3] - point[0:3])) == 1:
                return parallel.tolist().index(list(point))

    def _get_first_right_meridian_point_index(self, parallel, sorted_parallel, meridian_point):
        """
        This method assumes the use of spheres as iter_surface.
        """

        tan = cross(meridian_point[3:6],
                    meridian_point[0:3] / linalg.norm(meridian_point[0:3]))

        for point in sorted_parallel:
            if sign(dot(tan, meridian_point[0:3] - point[0:3])) == -1:
                return parallel.tolist().index(list(point))

    def compute_rays_from_grid(self, grid, grid_to_trajectories = None, points_to_num = None):
        """ This method compute the rays that belongs of the given grid. Note that you must create grid_to_trajectories
        first.

        Args:
            grid: (int) number of the grid
        """

        if grid_to_trajectories is None:
            grid_to_trajectories = self.load_grid_to_trajectories()
        parallels, borders = grid_to_trajectories[grid]
        rays = self.compute_rays_from_parallels(parallels, borders, points_to_num)
        return rays

    def compute_rays_from_parallels(self, parallels, borders = None, points_to_num = None):
        """ This method gets a list of point numbers (parallel db format) and
        create a list of rays (point-normal) with possible border points.
        It removes empty borders and empty parallels, but if both are empty,
        it will return an empty list.

        Args:
            parallels: list of [ list of (point numbers)]
            borders: If present, must be a list shape (N,2) with begin and end of each parallel
        """

        rays = []
        if points_to_num is None:
            points_to_num = self.load_points_to_num()

        ntp = self.get_sorted_points(points_to_num)
        parallels = copy.deepcopy(parallels)
        blade = self.load_blade()

        for i, parallel in enumerate(parallels):
            if len(parallel) == 0:
                rays += [[]]
                continue
            if self.verticalized == False:
                model = blade.select_model(ntp[parallel[0]])
                ray = map(lambda x: blade.compute_ray_from_point(ntp[x], model), parallel)
            else:
                ray = map(lambda x: blade.compute_ray_from_point(ntp[x]), parallel)

            if borders is not None:
                if len(borders[i][0]) > 0:
                    ray.insert(0,blade.compute_ray_from_point(borders[i][0]))

                if len(borders[i][1]) > 0:
                    ray.append(blade.compute_ray_from_point(borders[i][1]))
            rays.append(ray)
            
        return rays

    def _check_line(self, line, grid_bases, line_grid, line_grid_dist, threshold):
        """ Verify if given line (rail) can coat given grids.
        
        Args:
            line: rail configuration. ((x1,y1),(x2,y2))
            grid_bases: (dict) relates grid to bases (bases that can coat given grid).
            line_grid: (dict) relates lines (rails) to coatable grids.
            line_grid_dist: (dict) relates rail/grid to dist, a score for that line/grid combination.
            threshold: radius of each base.
        """

        x1 = line[0]
        x2 = line[1]
        grid_dist = dict()
        for grid, bases in grid_bases.iteritems():
            bases = list(bases)
            point_near, distance, distance_str = mathtools.distance_line_bases(
                x1, x2, bases, threshold)
            if distance != None:
                line_grid[line] = line_grid.get(line, set()) | set([grid])
                grid_dist[grid] = [point_near, distance, distance_str]
        line_grid_dist[line] = grid_dist
        return line_grid, line_grid_dist

    def compute_rail_configurations(self, lines, threshold):
        """ Compute rails configurations (lines) for all dbs. The method saves line_grid and line_grid_dist files.
            Returns:
                line_grid: (dict) relates lines (rails) to coatable grids. Call line_grid[line].
                line_grid_dist: (dict) some rail configuration (line) may not be the best to coat a specific grid.
                This dictionary relates rail/grid to dist, a score for that line/grid combination.
                Call line_grid_dist[line][grid].

        Args:
            lines: (list) list of tuples. All possible combinations of possible rail positions.
            [[((x1,y1),(x2,y2)),],[((x1,y1),(x2,y2)),((x3,y3),(x4,y4))]]
            threshold: radius of each base
        """

        dbs = self.info.findall('db')
        for dbi in dbs:
            line_grid = dict()
            line_grid_dist = dict()
            db_path = join(self.path, dbi.find('path').text)
            grid_bases = load_pickle(join(db_path, 'grid_bases.pkl'))
            for line in lines:
                line_grid, line_grid_dist = self._check_line(
                    line, grid_bases, line_grid, line_grid_dist, threshold)

            rail_path = join(db_path, 'rails')
            try:
                makedirs(rail_path)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise
            save_pickle(line_grid, join(rail_path, 'line_grid.pkl'))
            save_pickle(line_grid_dist, join(rail_path, 'line_grid_dist.pkl'))
        return

    def _select_db(self, grids_to_coat):
        """ Method returns combination of dbs to coat given grids. For example, if you want to coat grids [50,55] and
        you have 3 dbs: ['db_0','db_45','db_-45'], and only dbs ['db_0','db_45'] can coat both, the method will return
        [('db_0',),('db_45',)]. If only the combination ('db','db_45') can coat it (e.g. grid 50 is coatable by 'db_0'
        and grid 55 is coatable by 'db_45'), then the method will return: [('db_0','db_45')].

        Args:
            grids_to_coat: (list) grids to coat
        """

        db_grids = self.get_dbs_grids()

        for db, grids in db_grids.iteritems():
            if len(grids & set(grids_to_coat)) == 0:
                db_grids.pop('db',None)

        for ncombination in range(1, len(db_grids) + 1):
            feasible_combination = []
            for dbs in combinations(db_grids.keys(), ncombination):
                coatable = set([])
                for dbi in dbs:
                    coatable = coatable | db_grids[dbi]
                for grid in grids_to_coat:
                    if grid not in coatable:
                        break
                else:
                    feasible_combination.append(dbs)
            if len(feasible_combination) > 0:
                return feasible_combination
        return

    def _compute_minimal_lines(self, line_grid, grid_nums):
        lines = line_grid.keys()
        sol = []
        for i in range(1, len(lines)):
            for line_comb in combinations(lines, i):
                line_union = mathtools.union_line_grids(line_grid, line_comb)
                if len(grid_nums - line_union) == 0:
                    sol.append(line_comb)
            if len(sol) > 0:
                return sol
        return

    def _compute_psalpha_ordered_by_score(self, line_combs, grids_intersection, line_grid, line_grid_dist,
                                          criteria='sum'):
        """
        Given rail combinations, lines solutions to coat specific grids, 
        line_grid_dist dictionary, and criteria, it returns the ordered psalpha by score (criteria)
        of the rails.
        """

        min_str = []
        sum_str = []
        criteria = {'min': 1, 'sum': 0}[criteria]
        psalphas = dict()
        counter = 0
        for line_comb in line_combs:
            counter += 1
            distance_str_total = 0
            distance_str_min = 1000
            psalphas[line_comb] = dict()
            grids = grids_intersection
            for line in line_comb:
                grids_intersec = grids & line_grid[line]
                psalphas[line_comb][line] = dict()
                for grid in grids:
                    x1 = line[0][0]
                    y1 = line[0][1]
                    point_near, distance, distance_str = line_grid_dist[line][grid]
                    p = mathtools.closest_point_line_3d(array(line[0]), array(line[1]), point_near)
                    psalphas[line_comb][line][grid] = (
                    x1, sign(p[1]) * linalg.norm(p - line[0]), sign(p[1]) * atan2(x1 - p[0], abs(p[1] - y1)))
                    point_near, distance, distance_str = line_grid_dist[line][grid]
                    distance_str_total += distance_str
                    distance_str_min = min(distance_str_min, distance_str)
                grids = grids - grids_intersec
            sum_str.append(distance_str_total)
            min_str.append(distance_str_min)

        values = zeros((len(line_combs), 2))
        values[:, 0] = sum_str
        values[:, 1] = min_str
        ordered_sol = [x for (y, x) in sorted(zip(-values[:, criteria], line_combs))]
        return ordered_sol, psalphas

    def get_rail_configurations(self, grids_to_coat, criteria='sum'):
        """ This method returns rail solutions.

        Args:
            grids_to_coat: (list) grids to be coated.
            criteria: criteria to choose best placement
        """
        feasible_combinations = self._select_db(grids_to_coat)
        if feasible_combinations == None:
            print "Grids are not coatable"
            return

        db_grids = self.get_dbs_grids()
        sol_dict = dict()
        coated_grids = dict()
        psalphas = dict()
        for fcomb in feasible_combinations:
            coated_grids[fcomb] = dict()
            psalphas[fcomb] = dict()
            sol_dict[fcomb] = dict()
            grids = set(grids_to_coat)
            for dbi in fcomb:
                line_grid = load_pickle(join(dbi, 'rails', 'line_grid.pkl'))
                line_grid_dist = load_pickle(join(dbi, 'rails', 'line_grid_dist.pkl'))
                grids_intersection = grids & db_grids[dbi]
                line_combs = self._compute_minimal_lines(line_grid, grids_intersection)
                line_combs, psalpha = self._compute_psalpha_ordered_by_score(
                    line_combs, grids_intersection, line_grid, line_grid_dist, criteria)
                coated_grids[fcomb] = grids_intersection
                psalphas[fcomb][dbi] = psalpha
                sol_dict[fcomb][dbi] = line_combs
                grids = grids - db_grids[dbi]
        return sol_dict, psalphas, coated_grids, feasible_combinations

    def get_rail_configuration_n(self, grids_to_coat, criteria='sum', n=0):
        sol_dict, psalphas, coated_grids, feasible_combinations = self.get_rail_configurations(grids_to_coat, criteria)
        n_psas = []
        for fcomb in feasible_combinations:
            n_psa = []
            for dbi in fcomb:
                dbs = self.info.findall('db')
                for db in dbs:
                    if db.find('name').text == split(dbi)[1]:
                        T = self._extract_T(db)
                        T = atan2(-T[1,2],T[1,1])+pi/2
                try:
                    line_comb = sol_dict[fcomb][dbi][n]
                except IndexError:
                    break
                else:
                    for line in psalphas[fcomb][dbi][line_comb].keys():
                        for grid, psa in psalphas[fcomb][dbi][line_comb][line].iteritems():
                            rp = rail_place.RailPlace(psa)
                            xyz = rp.getTransform().flatten()
                            n_psa.append([list(psa),list(xyz)])
                break

            n_psa = [T, n_psa]
            n_psas.append(n_psa)
        return n_psas, feasible_combinations

    ################################################################################
    ############################## CLEAR METHODS ###################################
    ################################################################################

    def remove_point(self, points):
        """
        Remove given points from db.
        
        Keyword arguments:
        points (or rays) -- list of points to be removed
        """

        db = self.load_db()
        points_to_num = self.load_points_to_num()

        for point in points:
            try:
                key_num = points_to_num[tuple(round(point[0:3], 9))]
            except KeyError:
                continue

            db.pop(key_num, None)

        try:
            save_pickle(db, join(self.db_main_path, 'db.pkl'))
        except IOError:
            None

        return

    def clear_db(self):
        """
        Clear the db.
        """

        db = self.load_db()
        for key, value in db.iteritems():
            db[key] = set()
        try:
            save_pickle(db, join(self.db_main_path, 'db.pkl'))
        except IOError:
            raise
        return db

    def clear_visited_bases(self):
        """
        Clear the visited_bases.
        """

        visited_bases = self.load_visited_bases()
        for key in visited_bases:
            visited_bases[key] = False
        save_pickle(visited_bases, join(self.db_main_path, 'visited_bases.pkl'))
        return

    ################################################################################
    ############################### PLOT METHODS ###################################
    ################################################################################

    def plot_points_db(self, vis, scale=1):
        """
        Method to plot db points.

        keyword arguments:
        vis -- visualizer object.
        scale -- number_of_points/scale will be plotted.
        """

        db = self.load_db()
        index = map(int, random.uniform(0, len(db) - 1, int(len(db) * 1.0 / scale)))
        points_num = array(db.keys())[index]

        points_tuple = self.get_sorted_points()

        for i in points_num:
            _ = vis.plot(points_tuple[i], 'points_db')
        return

    def plot_bases_db(self, vis):
        """
        Method to plot db bases.

        keyword arguments:
        vis -- visualizer object.
        """

        db = self.load_db()
        db_bases = self.get_sorted_bases()
        bases = self.get_bases(db)
        for base in bases:
            rp = rail_place.RailPlace(db_bases[base])
            vis.plot(rp.getXYZ(self.turbine.config), 'base', (0, 0, 1))
        return

    def plot_grid(self, grid_num, vis):
        """ Method to plot grid.

        Args:
            grid_num: grid to be plotted.
            vis: visualizer object.
        """

        grid_to_trajectories = self.load_grid_to_trajectories()
        trajectories, borders = grid_to_trajectories[grid_num]
        ntp = self.get_sorted_points()
        for trajectory in trajectories:
            for point in trajectory:
                _ = vis.plot(ntp[point], 'p', color=(1, 0, 0))
        _ = vis.plot_lists(borders, 'p', color=(1, 0, 0))
        return

    def plot_grid_coat(self, vis, grid_num, base):
        """
        Method to plot coatable and non-coatable points of the grid.

        keyword arguments:
        grid_num -- grid to be plotted (int).
        vis -- visualizer object.
        base -- base (DB) for the robot (int).
        """

        non_coatable = []
        coatable = []
        grid_to_trajectories = self.load_grid_to_trajectories()
        trajectories, _ = grid_to_trajectories[grid_num]
        main_db = self.load_db()
        ntp = self.get_sorted_points()
        counter = 0
        for trajectory in trajectories:
            for point in trajectory:
                try:
                    if base not in main_db[point]:
                        non_coatable.append(ntp[point])
                    else:
                        coatable.append(ntp[point])
                except KeyError:
                    continue
        p = vis.plot(non_coatable, 'noncoat', color=(1, 0, 0))
        p = vis.plot(coatable, 'coat', color=(0, 0, 1))
        return non_coatable

    def plot_convex_grid(self, threshold, grid_num):
        """
        Method to plot convex-hull of bases that can coat a given grid.

        Keyword arguments:
        threshold -- minimum score to consider that the grid is coated.
        grid_num -- grid to coat.
        """

        import matplotlib.pyplot as plt
        bases, scores = self.base_grid_coating(grid_num)
        xy = []
        feasible_bases = []
        for i, base in enumerate(bases):
            if abs(scores[i]) >= threshold:
                rp = rail_place.RailPlace(base)
                xyz = rp.getXYZ(self.turbine.config)
                feasible_bases.append([scores[i], base])
                xy.append([xyz[0], xyz[1]])
            else:
                break

        if len(xy) > 0:
            fig = plt.figure()
            try:
                hull2D = ConvexHull(xy)
                plt = mathtools.plot_hull(xy, hull2D, plt)
            except:
                plt.scatter(array(xy)[:, 0], array(xy)[:, 1])
            fig.savefig(join(self.db_main_path, str(grid_num) + '.pdf'))
            plt.close()
        return
