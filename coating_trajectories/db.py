from path_filters import filter_trajectories, side_filter
import planning
from numpy import save, load, random, array, linspace, cross
from numpy import sign, dot, linalg, sum, zeros
from os.path import basename, splitext, join, exists, isfile
from os import makedirs, listdir
import copy
import rail_place
from colorsys import hls_to_rgb
import cPickle
import mathtools
import errno


class NoDBFound(Exception):    
    def __init__(self, value):
        Exception.__init__(self,"No DataBase named " + value + " found.")

class DB:
    """
    DB for robot's base position and coated points.

    Keyword arguments:
    path -- is the folder to store all data bases;
    blade -- BladeModeling object;
    """

    def __init__(self, path, blade=None):

        self.path = path
        self.create_db(blade)
        del blade
        
        if not exists(self.path):
            makedirs(self.path)
       
    def save_db_pickle(self, obj, name ):
        """
        Save a general file in pickle format.

        Keyword arguments:
        obj -- object to be saved.
        name -- name of the file ('example/file.pkl').
        """
        
        with open(name, 'wb') as f:
            cPickle.dump(obj, f, cPickle.HIGHEST_PROTOCOL)
        return

    def load_db_pickle(self, name):
        """
        Load a general pickle file.

        Keyword arguments:
        name -- name of the file ('example/file.pkl').
        """
        
        with open(name, 'rb') as f:
            return cPickle.load(f)

    def load_db(self):
        """
        Load main database num:num. The first num represents the points and
        the second represents the bases.
        """

        name = join(self.path,'fixed_db','db.pkl')
        return self.load_db_pickle(name)

    def create_db(self, blade):
        """
        Create the database.

        Keyword argument:
        blade -- the blade object. The database point should be created
        with the parallel points.
        """

        if not exists(join(self.path,'fixed_db')):
            makedirs(join(self.path,'fixed_db'))
        db = dict()
        db_points_to_num = dict()
        
        try:
            db = self.load_db()
        except IOError:
            if blade is None:
                db_points_to_num = dict()    
            else:
                counter = 0
                for trajectory in blade.trajectories:
                    for point in trajectory:
                        db[counter] = set()
                        counter+=1
            self.save_db_pickle(db, join(self.path,'fixed_db','db.pkl'))
        del db
            
        try:
            db_points_to_num = self.load_db_points_to_num()
        except IOError:
            if blade is None:
                db_points_to_num = dict()
            else:
                counter = 0
                for trajectory in blade.trajectories:
                    for point in trajectory:
                        db_points_to_num[tuple(point[0:3])] = counter
                        counter+=1
            self.save_db_pickle(db_points_to_num, join(self.path,'fixed_db','db_points_to_num.pkl'))
        del db_points_to_num

        try:
            db_bases_to_num = self.load_db_bases_to_num()
        except IOError as error:
            if error.errno == errno.ENOENT:
                raise NoDBFound('db_bases_to_num.pkl')
            else:
                raise
        del db_bases_to_num
           
        return

    def load_db_points_to_num(self):
        """
        Load points_to_num database points:num. Points are tuples (x,y,z) and
        the num maps the points (counter to reduce database complexity).
        """

        path = join(self.path, 'fixed_db', 'db_points_to_num.pkl')
        return self.load_db_pickle(path)

    def load_db_grid_to_mp(self):
        """
        Load grid_to_mp database num:[(m1,m2),(p1,p2)].
        Grids are numbers and the mp are
        meridian/parallels index.
        """

        path = join(self.path, 'fixed_db', 'db_grid_to_mp.pkl')
        return self.load_db_pickle(path)

    def load_db_grid_to_bases(self):
        """
        Load grid_to_bases database num:[bases].
        Grids are numbers and the bases are lists of PSAlpha tuples.
        """

        path = join(self.path, 'fixed_db', 'db_grid_to_bases.pkl')
        return self.load_db_pickle(path)

    def load_db_grid_to_trajectories(self):
        """
        Load grid_to_trajectories database num:[trajectories, borders].
        Grids are numbers and the bases are lists of trajectories (num, db format)
        and borders (tuples Nx3).
        """

        path = join(self.path, 'fixed_db', 'db_grid_to_trajectories.pkl')
        return self.load_db_pickle(path)

    def load_db_bases_to_num(self):
        """
        Load bases_to_num database bases:num. Bases are tuples railplace and
        the num maps the bases (counter to reduce database complexity).
        """

        if not exists(join(self.path,'fixed_db')):
            makedirs(join(self.path,'fixed_db'))

        path = join(self.path,'fixed_db','db_bases_to_num.pkl')
        
        return self.load_db_pickle(path)

    def load_db_visited_bases(self):
        """
        Load visited_bases database bases:bool. Bases are tuples PSAlpha and
        the bool show if it was already computed.
        """

        path = join(self.path,'fixed_db','db_visited_bases.pkl')
        return self.load_db_pickle(path)

    def load_grid_meridian(self):
        """
        Load grid_meridians. The meridians used to make the grid.
        """

        path = join(self.path, 'fixed_db', 'meridians.pkl')
        return self.load_db_pickle(path)

    def load_grid_parallel(self):
        """
        Load grid_parallel. The parallels used to make the grid.
        """

        path = join(self.path, 'fixed_db', 'parallels.pkl')
        return self.load_db_pickle(path)

    def get_sorted_bases(self):
        """
        Return the sorted bases -- tuples PSAlpha.
        """

        btn = self.load_db_bases_to_num()
        return [ b for (v,b) in sorted(zip(btn.values(),btn.keys()))]

    def get_sorted_points(self):
        """
        Return the sorted points -- tuples (x,y,z).
        """

        ptn = self.load_db_points_to_num()
        return [ b for (v,b) in sorted(zip(ptn.values(),ptn.keys()))]

    def get_bases(self, db):
        bases = set()
        for val in db.values():
            bases = val|bases
        return bases

    def merge_db(self, db_file, main_db):
        """
        Merge new db_file to the main database.

        keyword arguments:
        db_file -- file to be merged;
        main_db -- main database file;
        """

        for key, value in db_file.iteritems():
            if type(value) == int:
                value = set([value])
            main_db[key] = main_db.get(key,set()) | value
        return main_db
        
    def merge_db_directory(self, path):
        """
        Method to merge all databases in a specific folder with the
        main database.

        keyword arguments:
        path -- where the databases are
        """

        db = self.load_db()
        onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
        N = len(onlyfiles)
        index = 0
        for afile in onlyfiles:
            filename, file_extension = splitext(afile)
            if file_extension == '.pkl':
                index+=1
                print('iter %3i\tfile = %s' % (index, filename)) 
                db = self.merge_db(self.load_db_pickle(join(path,afile)), db)
        self.save_db_pickle(db, join(self.path,'fixed_db','db.pkl'))
        return

    def generate_db(self, turbine, blade, rail_position, minimal_number_of_points_per_trajectory = 2): 
        """
        Function that store feasible trajectories for a given robot's base position.
        The robot must be placed before this method call.
        Each feasible point in a trajectory is mapped to a PSAlpha that can coat it.
        Therefore, the DB returns
        a dictionary where points are keys and rail positions are values, e.g.:
             {     (1.0, 1.1, 1.3):                   (2.4   ,  1.2  ), ( 2.4  ,  1.3  )}
             {point( x ,  y ,  z ): rail poisitions   (rail_1, rail_2), (rail_1, rail_2)}

        Keyword arguments:
        turbine -- turbine object
        rail_position -- PSAlpha, base_number or rail_place object. 
        minimal_number_of_points_per_trajectory -- trajectories smaller then
        this number of points are croped
        """

        if isinstance(rail_position, rail_place.RailPlace):
            rail_position = rail_place.getPSAlpha()
        if type(rail_position) == tuple:
            db_bases_to_num = self.load_db_bases_to_num()
            base = db_bases_to_num[rail_position]
            del db_bases_to_num
        elif type(rail_position) == int:
            base = rail_position
        else: raise TypeError('rail_position is not valid')
	
        db_points_to_num = self.load_db_points_to_num()
        db = dict()

        filtered_trajectories = filter_trajectories(turbine, blade.trajectories,2)
        for filtered_trajectory in filtered_trajectories:
            for filtered_trajectory_part in filtered_trajectory:
                evaluated_points = 0
                while evaluated_points < len(filtered_trajectory_part):
                    try:
                        lower, _, _ = planning.compute_first_feasible_point(
                            turbine,
                            filtered_trajectory_part[evaluated_points:],
                            blade.trajectory_iter_surface)
                        evaluated_points = evaluated_points + lower
                    except ValueError: 
                        evaluated_points = len(filtered_trajectory_part)
                        continue

                    joint_solutions = planning.compute_robot_joints(turbine,
                                                                    filtered_trajectory_part,
                                                                    evaluated_points,
                                                                    blade.trajectory_iter_surface)
                    upper = evaluated_points + len(joint_solutions)
                    if upper-evaluated_points > minimal_number_of_points_per_trajectory:
                        for point in filtered_trajectory_part[evaluated_points:upper]:
                            db[db_points_to_num[tuple(point[0:3])]] = set([base])
                    evaluated_points = upper+1
        return db

    def invert_db(self, db, parallels, regions):
    # Regions as list of [list of tuple (parallel_index, begin_index, end_index)]
        psa_db = list()
        
        for region in regions:
            inverse_set = set()
            
            for parallel_index, begin_index, end_index in region:
                for point in parallels[parallel_index, begin_index:end_index]:
                    inverse_set |= db[point]

            psa_db.append(inverse_set)
        return

    def get_bases_trajectories(self, trajectories):
        """
        Method returns the robot bases (tuples PSAlpha)
        that can coat all set of points in trajectories and db.
        If a point in trajectories does not belong to the db, it is
        skipped.

        keyword arguments:
        trajctories -- points_num (db format) to be coated
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

    def get_best_bases_trajectories(self, trajectories):
        """
        Method returns the bases that can coat not all, but almost
        all of points in side the grid. This is important to know
        because some points should not be coated (as points on top of
        the blade)

        keyword arguments:
        trajctories -- points_num (db format) to be coated
        """

        db = self.load_db()
        bases_tuple = self.get_sorted_bases()
        score = zeros(len(bases_tuple))
        N = 0
        
        for trajectory in trajectories:
            N+=len(trajectory)
            for point in trajectory:
                try:
                    score[list(db[point])]-=1
                except KeyError:
                    pass
        best_bases = [x for (y,x) in sorted(zip(score,bases_tuple))]
        return best_bases, -array(sorted(score))*1.0/N
        
    def make_grid(self, blade, number_of_meridians = 12, number_of_parallels = 6, init_parallel = 17 ):
        """
        Make a grid in the blade with parallels and meridians.
        Meridians and parallels are evenly spaced.
        The parallel used to compute the meridians is in the middle.
        
        keyword arguments:
        blade -- blade object.
        number_of_meridians -- number of meridians to be computed.
        number_of_parallels -- number of parallels.
        init_parallel -- initial parallel. This argument is
        important to remove the lip from the grid.
        """
        
        parallel = blade.trajectories[int(len(blade.trajectories)/2)]
        meridians = blade.draw_meridians(parallel, 1e-3, number_of_meridians)
        
        parallels = []
        list_index = linspace(init_parallel,len(blade.trajectories),number_of_parallels).astype(int)
        list_index[-1] -= 2
        for i in list_index:
            parallels.append(blade.trajectories[i])
        return meridians, parallels

    def create_db_grid(self, blade):
        """
        Autonomously create the db_grid after the make_grid process.
        The method will get adjacent meridians and parallels, building
        the grids. Some of the grids may not make sense, thus the user
        should analyse, add and remove grids.
        
        keyword arguments:
        blade -- blade object.
        """
        
        db_grid_to_mp = dict()
        db_grid_to_bases = dict()
        db_grid_to_trajectories = dict()
        meridians = self.load_grid_meridian()
        parallels = self.load_grid_parallel()

        counter = 0
        for i in range(0,len(meridians)):
            for j in range(0,len(parallels)-1):
                grid = [(i, (i+1)%len(meridians)), (j,j+1)]
                print('iter %3i ' % (counter))
                trajectories_in_grid, borders = self.get_points_in_grid(
                    blade, [meridians[grid[0][0]], meridians[grid[0][1]]],
                    [parallels[grid[1][0]], parallels[grid[1][1]]])
                bases = self.get_bases_trajectories(trajectories_in_grid)

                db_grid_to_mp[counter] = grid
                db_grid_to_bases[counter] = bases 
                db_grid_to_trajectories[counter] = [trajectories_in_grid, borders]
                counter+=1

        try:
            self.save_db_pickle(db_grid_to_mp, join(self.path,'fixed_db','db_grid_to_mp.pkl'))
        except IOError: None

        try:
            self.save_db_pickle(db_grid_to_bases, join(self.path,'fixed_db','db_grid_to_bases.pkl'))
        except IOError: None

        try:
            self.save_db_pickle(db_grid_to_trajectories, join(self.path,'fixed_db','db_grid_to_trajectories.pkl'))
        except IOError: None
       
        return

    def create_db_grid_to_bases(self, db_grid_to_trajectories=None):
        """
        Compute bases to coat the grids
        
        keyword arguments:
        blade -- blade object.
        T -- homogenous transform matrix to rotate the blade.
        """

        if db_grid_to_trajectories is None:
            db_grid_to_trajectories = self.load_db_grid_to_trajectories()
        db_grid_to_bases = dict()

        for key, value in db_grid_to_trajectories.iteritems():
            trajectories, borders = value
            bases = self.get_bases_trajectories(trajectories)
            db_grid_to_bases[key] = bases 

        try:
            self.save_db_pickle(db_grid_to_bases, join(self.path,'fixed_db','db_grid_to_bases.pkl'))
        except IOError: None
        return
                           

    def get_points_in_grid(self, blade, meridian, parallel):
        """
        Get parallels that belong to a grid, between given meridians and
        given parallels.
        
        keyword arguments:
        blade -- blade object.
        meridian -- tuple 1x2, First meridian must be on the right w.r.t. the second meridian.
        It means that points[1] (y) of points meridian1 > points[1] (y) of points meridian2.  
        parallel -- tuple 1x2. Not ordered.
        """

        meridian1, meridian2 = meridian[0], meridian[1]
        parallel1, parallel2 = parallel[0], parallel[1]
        db_points_to_num = self.load_db_points_to_num()
        
        parallel_index_1 = int((blade.trajectory_iter_surface._Rn0 -
                                linalg.norm(parallel1[0][0:3]))/
                               blade.trajectory_iter_surface.coatingstep)
        parallel_index_2 = int((blade.trajectory_iter_surface._Rn0 -
                                linalg.norm(parallel2[0][0:3]))/
                               blade.trajectory_iter_surface.coatingstep)
        init = min(parallel_index_1,parallel_index_2)
        end = max(parallel_index_1,parallel_index_2)+1

        trajectories_in_grid = []

        def get_point_value(point):
            try:
                return db_points_to_num[point]
            except KeyError:
                return None

        borders = []
        for i in range(init,end):
            trajectory_in_grid = []
            parallel = array(blade.trajectories[i])
            blade.trajectory_iter_surface.find_iter(parallel[0])
            p1, sorted_parallel1 = self._closest_meridian_point(meridian1, parallel, blade)
            p2, sorted_parallel2 = self._closest_meridian_point(meridian2, parallel, blade)
            parallel1 = mathtools.direction_in_halfplane(parallel,p1[3:6])
            parallel2 = mathtools.direction_in_halfplane(parallel,p2[3:6])
            p1, sorted_parallel1 = self._closest_meridian_point(meridian1, parallel1, blade)
            p2, sorted_parallel2 = self._closest_meridian_point(meridian2, parallel2, blade)
            index_left = self._get_first_left_meridian_point_index(parallel, sorted_parallel1, p1)
            index_right = self._get_first_right_meridian_point_index(parallel, sorted_parallel2, p2)
            
            if abs(index_right - index_left)%(len(parallel)-1) == 1:
                pass
            elif index_left <= index_right:
                trajectory_in_grid += filter(lambda x: x is not None,
                                             map(get_point_value,
                                                 map(tuple, array(parallel)[index_left:index_right+1, 0:3].tolist())))
            else:
                trajectory_in_grid += filter(lambda x: x is not None,
                                             map(get_point_value,
                                                 map(tuple, array(parallel)[index_left:, 0:3].tolist())))+\
                                      filter(lambda x: x is not None,
                                             map(get_point_value,
                                                 map(tuple, array(parallel)[:index_right+1, 0:3].tolist())))   
            trajectories_in_grid.append(trajectory_in_grid)
            borders.append([tuple(p1[0:3]),tuple(p2[0:3])])
        return trajectories_in_grid, borders

    def _closest_meridian_point(self, meridian, parallel, blade):
        min_dist = 100
        closest_meridian_point = []
        sorted_parallel = []
        dist_list = []
        for meridian_point in meridian:
            dist = sum((parallel[:,0:3]-meridian_point[0:3])*
                       (parallel[:,0:3]-meridian_point[0:3]),1)
            if min(dist) <= min_dist:
                closest_meridian_point = meridian_point
                min_dist = min(dist)
                dist_list = dist
        sorted_parallel = [x for (y,x) in sorted(zip(dist_list,parallel))]
        model = blade.select_model(closest_meridian_point)
        return mathtools.curvepoint(model,blade.trajectory_iter_surface,closest_meridian_point[0:3]), sorted_parallel

    def _get_first_left_meridian_point_index(self, parallel, sorted_parallel, meridian_point):
        """
        This method assumes the use of spheres as iter_surface.
        """
        
        tan = cross(meridian_point[3:6],
                    meridian_point[0:3]/linalg.norm(meridian_point[0:3]))
        for point in sorted_parallel:
            if sign(dot(tan,meridian_point[0:3]-point[0:3])) == 1:
                return parallel.tolist().index(list(point))

    def _get_first_right_meridian_point_index(self, parallel, sorted_parallel, meridian_point):
        """
        This method assumes the use of spheres as iter_surface.
        """
        
        tan = cross(meridian_point[3:6],
                    meridian_point[0:3]/linalg.norm(meridian_point[0:3]))

        for point in sorted_parallel:
            if sign(dot(tan,meridian_point[0:3]-point[0:3])) == -1:
                return parallel.tolist().index(list(point))

    def compute_rays_from_parallels(self, blade, parallels, borders = None):
        """
        This method gets a list of point numbers (parallel db format) and
        create a list of rays (point-normal) with possible border points.
        It removes empty borders and empty parallels, but if both are empty,
        it will return an empty list.

        Keyword arguments:
        blade -- a blade_modeling object.
        parallels -- list of [ list of (point numbers)]
        borders -- If present, must be a list shape (N,2) with begin and end of each parallel
        """

        rays = []
        ntp = self.get_sorted_points()

        def get_ray(model,point):
            df = model.df(point)
            df = df/linalg.norm(df)
            return array(list(point)+list(df))
            
        if borders is None:
            for parallel in parallels:
                if len(parallel)==0:
                    rays+=[[]]
                    continue
                model = blade.select_model(ntp[parallel[0]])
                rays += [ map( lambda x: get_ray(model,ntp[x]), parallel ) ]

        else:
            for i in range(len(parallels)):
                traj = []
                if len(borders[i][0])>0:
                    model = blade.select_model(borders[i][0])
                    traj.append(get_ray(model,borders[i][0]))

                if len(parallels[i])>0:
                    model = blade.select_model(ntp[parallels[i][0]])
                    traj += map( lambda x: get_ray(model,ntp[x]), parallels[i])
                    
                if len(borders[i][1])>0:
                    model = blade.select_model(borders[i][1])
                    traj.append(get_ray(model,borders[i][1]))
                    
                rays.append(traj)
        return rays
                
    def bases_validation(self, parallels, bases, turbine, blade):
        """
        Validate bases for given parallels.
        
        Keyword arguments:
        parallels -- list of rays (N,6)
        bases -- [PSAlpha] (N,3), list of tuples
        """
        
        feasible_bases = []
        for base in bases:
            rp = rail_place.RailPlace(base)
            turbine.place_rail(rp)
            turbine.place_robot(rp)
            for parallel in parallels:
                joint_solutions = planning.compute_robot_joints(
                    turbine, parallel, 0, blade.trajectory_iter_surface)
                if len(joint_solutions) != len(parallel):
                    break
            else:
                feasible_bases.append(base)
            continue
        return feasible_bases

    def remove_point(self, points):
        """
        Remove given points from db, db_points_to_num and
        db_grid_to_trajectories.
        
        Keyword arguments:
        points (or rays) -- list of points to be removed
        """

        db_grid_to_trajectories = self.load_db_grid_to_trajectories()
        db = self.load_db()
        db_points_to_num = self.load_db_points_to_num()
        
        for point in points:
            key_point = tuple(point[0:3])     
            try:
                key_num = db_points_to_num[key_point]
            except KeyError:
                continue

            db.pop(key_num,None)
            for key, value in db_grid_to_trajectories.iteritems():
                trajectories, border = value
                for trajectory in trajectories:
                    if key_num in trajectory:
                        trajectory.remove(key_num)

        try:
            self.save_db_pickle(db_grid_to_trajectories, join(self.path,'fixed_db','db_grid_to_trajectories.pkl'))
        except IOError: None
        try:
            self.save_db_pickle(db_points_to_num, join(self.path,'fixed_db','db_points_to_num.pkl'))
        except IOError: None
        try:
            self.save_db_pickle(db, join(self.path,'fixed_db','db.pkl'))
        except IOError: None

        return
                
    def clear_db(self):
        """
        Clear the db.
        """

        db = self.load_db()
        for key, value in db.iteritems():
            db[key] = set()
        return db

    def clear_db_visited_bases(self):
        """
        Clear the db_visited_bases.
        """

        db = self.load_db_visited_bases()
        for key, value in db.iteritems():
            db[key] = False
        return db

    def plot_points_gradient(self, vis, scale = 1):
        """
        Method to plot points in a color gradient way regarding reachability by the bases.
        Red -> point is not reachable by any base
        Green -> point is well reachable.

        keyword arguments:
        vis -- visualizer object.
        scale -- number_of_points/scale will be plotted.
        """
        
        N = 0
        db = self.load_db()
        for key, value in db.iteritems():
            N = max(N,len(value))

        points_tuple = self.get_sorted_points()
        
        index = map(int,random.uniform(0,len(db)-1,int(len(db)*1.0/scale)))
        points_num = array(db.keys())[index]
        bases_num = array(db.values())[index]

        for i in range(0,len(points_num)):
            vis.plot(points_tuple[points_num[i]], 'points_gradient',
                     color = hls_to_rgb(len(bases_num[i])*1.0/(3*N),0.5,1))
        return

    def plot_points_db(self, vis, scale=1):
        """
        Method to plot db points.

        keyword arguments:
        vis -- visualizer object.
        scale -- number_of_points/scale will be plotted.
        """

        db = self.load_db()
        index = map(int,random.uniform(0,len(db)-1,int(len(db)*1.0/scale)))
        points_num = array(db.keys())[index]

        points_tuple = self.get_sorted_points()
        
        for i in points_num:
            vis.plot(points_tuple[i], 'points_db')
        return

    def plot_bases_db(self, vis, turbine):
        """
        Method to plot db bases.

        keyword arguments:
        vis -- visualizer object.
        scale -- number_of_points/scale will be plotted.
        """

        db = self.load_db()
        db_bases = self.get_sorted_bases()
        bases = self.get_bases(db)
        for base in bases:
            rp = rail_place.RailPlace(db_bases[base])
            vis.plot(rp.getXYZ(turbine.config),'base',(0,0,1))
        return

    def plot_grid(self, blade, grid_num, vis):
        db_grid_to_trajectories = self.load_db_grid_to_trajectories()
        trajectories, borders = db_grid_to_trajectories[grid_num]
        rays = self.compute_rays_from_parallels(blade, trajectories, borders)
        vis.plot_lists(rays, 'rays', color=(1,0,0))
        return
        
