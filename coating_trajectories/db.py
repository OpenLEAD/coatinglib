from path_filters import filter_trajectories, side_filter
import planning
from numpy import save, load, random, array, linspace, cross
from numpy import sign, dot, linalg, sum, zeros, round
from os.path import basename, splitext, join, exists, isfile, split, realpath
from os import makedirs, listdir
import copy
import rail_place
import cPickle
import mathtools
import errno
from lxml import etree as ET
from openravepy import matrixFromAxisAngle
import blade_modeling


class NoDBFound(Exception):    
    def __init__(self, value):
        Exception.__init__(self,"No DB file named " + value + " found.")

def save_pickle(obj, name ):
    """
    Save a general file in pickle format.

    Keyword arguments:
    obj -- object to be saved.
    name -- name of the file ('example/file.pkl').
    """
    
    with open(name, 'wb') as f:
        cPickle.dump(obj, f, cPickle.HIGHEST_PROTOCOL)
    return

def load_pickle(filename):
    with open(filename,'rb') as f:
        return cPickle.load(f)

def get_sorted(filename):
    b2n = load_db_pickle(filename)
    return [ b for (v,b) in sorted(zip(b2n.values(),b2n.keys()))]
    

class DB:
    """
    DB for robot's base position and coated points.

    Keyword arguments:
    path -- is the folder to store all data bases;
    blade -- BladeModeling object;
    """

    def __init__(self, path, turbine, db_name = None):

        self.path = path
        self.db_main_path = None
        self.T = matrixFromAxisAngle([0,0,0])
        self.info = None
        self.turb = turbine

        if not exists(self.path):
            makedirs(self.path)

        try:
            self.info = ET.parse(open(join(self.path,'info.xml')))
        except IOError as error:
            if error.errno == errno.ENOENT:
                raise NoDBFound('info.xml')

        if db_name!=None:
            dbs = self.info.findall('db')
            db_main = None
            for db in dbs:
                if db.find('name').text == db_name:
                    self.db_main_path = join(self.path,db.find('path').text)
                    T = db.find('transform').text
                    T = T.replace('\n',' ')
                    T = T.replace('\t',' ')
                    T = T.split(' ')
                    T[:] = [x for x in T if x != '']
                    m = []; mi = []
                    for t in T:
                        try:
                            mi.append(float(t))
                        except: None
                        if len(mi)==4:
                            m.append(mi)
                            mi=[]
                    m = array(m)
                    if m.shape==(4,4):
                        self.T = m
                    else:
                        raise SyntaxError("Invalid transform shape in info.xml")
                    break
            else:
                raise NoDBFound(db_name)
                
                
    def load_db(self):
        """
        Load main database num:num. The first num represents the points and
        the second represents the bases.
        """

        name = join(self.db_main_path,'db.pkl')
        return load_pickle(name)

    def load_points_to_num(self):
        """
        Load points_to_num database points:num. Points are tuples (x,y,z) and
        the num maps the points (counter to reduce database complexity).
        """

        path = join(self.path, 'points_to_num.pkl')
        ptn = load_pickle(path)
        real_ptn = dict()
        for key, value in ptn.iteritems():
            real_ptn[tuple(round(dot(self.T,[key[0],key[1],key[2],1])[0:3],9))] = value
        return real_ptn


    def load_grid_to_mp(self):
        """
        Load grid_to_mp database num:[(m1,m2),(p1,p2)].
        Grids are numbers and the mp are
        meridian/parallels index.
        """

        path = join(self.path, 'grid_to_mp.pkl')
        return load_pickle(path)

    def load_grid_to_trajectories(self):
        """
        Load grid_to_trajectories database num:[trajectories, borders].
        Grids are numbers and the bases are lists of trajectories (num, db format)
        and borders (tuples Nx3).
        """

        path = join(self.path, 'grid_to_trajectories.pkl')
        grid_to_trajectories = load_pickle(path)
        new_grid_to_trajectories = dict()
        for grid, value in grid_to_trajectories.iteritems():
            trajectories, borders = value
            new_grid_to_trajectories[grid] = [trajectories,mathtools.rotate_points(borders, self.T)]
            
        return new_grid_to_trajectories

    def load_bases_to_num(self):
        """
        Load bases_to_num database bases:num. Bases are tuples railplace and
        the num maps the bases (counter to reduce database complexity).
        """

        path = join(self.path,'bases_to_num.pkl')
        
        return load_pickle(path)

    def load_visited_bases(self):
        """
        Load visited_bases database bases:bool. Bases are tuples PSAlpha and
        the bool show if it was already computed.
        """

        path = join(self.db_main_path,'visited_bases.pkl')
        return load_pickle(path)

    def load_grid_meridian(self):
        """
        Load grid_meridians. The meridians used to make the grid.
        """

        path = join(self.path, 'meridians.pkl')
        return mathtools.rotate_trajectories(load_pickle(path),self.T)

    def load_grid_parallel(self):
        """
        Load grid_parallel. The parallels used to make the grid.
        """

        path = join(self.path, 'parallels.pkl')
        return mathtools.rotate_trajectories(load_pickle(path),self.T)

    def load_blade(self):
        """
        Function to load blade model.
        """

        folder = self.info.find('blade_model_filtered')
        if folder == None:
            folder = self.info.find('blade_model')
        if folder == None:
            raise NameError("No blade model found in info.xml")

        folder = folder.text
        xml_trajectories_path = join(self.path,folder,"trajectory/trajectory.xml")
        blade = blade_modeling.BladeModeling(self.turb, self.turb.blades[0])
        blade.load_trajectory(xml_trajectories_path)
        blade.trajectories = mathtools.rotate_trajectories(blade.trajectories,self.T)
        blade.rotate_models(self.T)
        return blade

    def load_blade_full(self):
        """
        Function to load blade model.
        """

        folder = self.info.find('blade_model')
        if folder == None:
            raise NameError("No blade model found in info.xml")

        folder = folder.text
        xml_trajectories_path = join(self.path,folder,"trajectory/trajectory.xml")
        blade = blade_modeling.BladeModeling(self.turb, self.turb.blades[0])
        blade.load_trajectory(xml_trajectories_path)
        blade.trajectories = mathtools.rotate_trajectories(blade.trajectories,self.T)
        blade.rotate_models(self.T)
        return blade

    def create_points_to_num(self):
        blade = self.load_blade()
        points_to_num = dict()

        counter = 0
        for trajectory in blade.trajectories:
            for point in trajectory:
                points_to_num[tuple(round(point[0:3],9))] = counter
                counter+=1
        save_pickle(points_to_num, join(self.path, 'points_to_num.pkl'))
        return 

    def create_db(self):
        """
        Create the database.

        Keyword argument:
        blade -- the blade object. The database point should be created
        with the parallel points.
        """

        try:
            points_to_num = self.load_points_to_num()
        except IOError:
            self.create_points_to_num()
        
        try:
            db = self.load_db()
        except IOError:
            db = dict()
            blade = self.load_blade()
            counter = 0
            for trajectory in blade.trajectories:
                for point in trajectory:
                    db[counter] = set()
                    counter+=1
            save_pickle(db, join(self.db_main_path,'db.pkl'))
        del db

        try:
            bases_to_num = self.load_bases_to_num()
        except IOError as error:
            if error.errno == errno.ENOENT:
                raise NoDBFound('bases_to_num.pkl')
            else:
                raise        
        return

    def get_sorted_bases(self):
        """
        Return the sorted bases -- tuples PSAlpha.
        """

        btn = self.load_bases_to_num()
        return [ b for (v,b) in sorted(zip(btn.values(),btn.keys()))]

    def get_sorted_points(self):
        """
        Return the sorted points -- tuples (x,y,z).
        """

        ptn = self.load_points_to_num()
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

    def merge_seg(self, seg, main_db):
        """
        Merge new seg to the main database.

        keyword arguments:
        seg -- segments to be merged;
        main_db -- main database file;
        """

        base = seg.keys()[0]
        for segments in seg[base]:
            for segment in segments:
                for point in segment:
                    main_db[point] = main_db.get(point,set()) | set([base])
        return main_db

    def create_db_from_segments(self, path, merge = 0):
        """
        Method to create db (num to num) with segments in given path.

        keyword arguments:
        path -- where the segments are.
        """

        if merge==0:
            db = dict()
            ptn = self.load_points_to_num()
            for key, value in ptn.iteritems():
                db[value] = set()
        else: db = self.load_db()
                  
        onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
        N = len(onlyfiles)
        index = 0
        for afile in onlyfiles:
            filename, file_extension = splitext(afile)
            if file_extension == '.pkl':
                print('file = %s' % (filename))
                try: 
                    db_file = load_pickle(join(path,afile))
                except EOFError:
                    continue
                db = self.merge_seg(db_file, db)        
        return db

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
        Method returns the bases in score order that can coat the given
        trajectories. This order is more important to know than absolut coating
        because some points may never be coated (as points on the top of
        the blade)

        keyword arguments:
        trajectories -- points_num (db format) to be coated
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
        
    def make_grid(self, number_of_meridians, number_of_parallels, init_parallel):
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

        blade = load_blade_full(self)
        parallel = blade.trajectories[int(len(blade.trajectories)/2)]
        meridians = blade.draw_meridians(parallel, 1e-3, number_of_meridians)
        
        parallels = []
        list_index = linspace(init_parallel,len(blade.trajectories),number_of_parallels).astype(int)
        list_index[-1] -= 1
        for i in list_index:
            parallels.append(blade.trajectories[i])
        return meridians, parallels

    def create_db_grid(self):
        """
        Autonomously create the db_grid after the make_grid process.
        The method will get adjacent meridians and parallels, building
        the grids. Some of the grids may not make sense, thus the user
        should analyse, add and remove grids.
        
        keyword arguments:
        blade -- blade object.
        """

        blade = load_blade_full()
        grid_to_mp = dict()
        db_grid_to_bases = dict()
        grid_to_trajectories = dict()
        meridians = self.load_grid_meridian()
        parallels = self.load_grid_parallel()

        counter = 0
        for i in range(0,len(meridians)):
            for j in range(0,len(parallels)-1):
                grid = [(i, (i+1)%len(meridians)), (j,j+1)]
                trajectories_in_grid, borders = self.get_points_in_grid(
                    blade, [meridians[grid[0][0]], meridians[grid[0][1]]],
                    [parallels[grid[1][0]], parallels[grid[1][1]]])
                bases = self.get_bases_trajectories(trajectories_in_grid)

                grid_to_mp[counter] = grid
                db_grid_to_bases[counter] = bases 
                grid_to_trajectories[counter] = [trajectories_in_grid, borders]
                counter+=1

        try:
            save_pickle(grid_to_mp, join(self.path,'grid_to_mp.pkl'))
        except IOError: None

        try:
            save_pickle(db_grid_to_bases, join(self.path,'db_grid_to_bases.pkl'))
        except IOError: None

        try:
            save_pickle(grid_to_trajectories, join(self.path,'grid_to_trajectories.pkl'))
        except IOError: None
       
        return

    def generate_db_joints(self, base_num, minimal_number_of_points_per_trajectory=3):
        """
        Compute joints for coating HARD COMPUTATION time
        outputs:
        db_base_to_joints = { base: { point_num: joints } }
        db_base_to_segs = { base: [trajectory1, trajectory2, ...] }
        trajectory = [part1, part2,...]
        part = [point1_num, point2_num, ...]
        
        """

        ntb = self.get_sorted_bases()
        base = ntb[base_num]
        rp = rail_place.RailPlace(base)
        self.turb.place_rail(rp)
        self.turb.place_robot(rp)
        if self.turb.check_rail_collision():
            return None, None
        if self.turb.check_robotbase_collision():
            return None, None

        ptn = self.load_points_to_num()
        blade = self.load_blade()
        db_base_to_joints = dict()  
        db_base_to_segs = dict()
        db_base_to_joints[base_num] = dict()
        db_base_to_segs[base_num] = list()

        # iterate each (filtered) parallel
        filtered_trajectories = filter_trajectories(self.turb, blade.trajectories, minimal_number_of_points_per_trajectory)
        for filtered_trajectory in filtered_trajectories:
            db_base_to_segs[base_num] += [[]]
            # iterate each part of trajectory
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
                    
                    if len(joint_solutions) >= minimal_number_of_points_per_trajectory:
                        db_base_to_segs[base_num][-1] += [[]]
                        # save point_num in each dictionary 
                        for point, joints in zip(filtered_trajectory_part[evaluated_points:upper], joint_solutions):
                            point_num = ptn[tuple(round(point[0:3],9))]
                            db_base_to_joints[base_num][point_num] = joints
                            db_base_to_segs[base_num][-1][-1] += [point_num]

                    # restart at end point
                    evaluated_points = upper
        return db_base_to_joints, db_base_to_segs                           

    def get_points_in_grid(self, meridian, parallel):
        """
        Get parallels that belong to a grid, between given meridians and
        given parallels.
        
        keyword arguments:
        blade -- blade object.
        meridian -- tuple 1x2, First meridian must be on the right w.r.t. the second meridian.
        It means that points[1] (y) of points meridian1 > points[1] (y) of points meridian2.  
        parallel -- tuple 1x2. Not ordered.
        """

        blade = self.load_blade()
        meridian1, meridian2 = meridian[0], meridian[1]
        parallel1, parallel2 = parallel[0], parallel[1]
        points_to_num = self.load_points_to_num()

        parallel_index_1 = 0
        parallel_index_2 = 0
        for i in range(0,len(blade.trajectories)):
            traj_list = [list(a) for a in blade.trajectories[i]]
            if list(parallel1[0]) in traj_list:
                parallel_index_1 = i
                break
        for i in range(0,len(blade.trajectories)):
            traj_list = [list(a) for a in blade.trajectories[i]]
            if list(parallel2[0]) in traj_list:
                parallel_index_2 = i
                break

        init = min(parallel_index_1,parallel_index_2)
        end = max(parallel_index_1,parallel_index_2)+1
        trajectories_in_grid = []

        def get_point_value(point):
            try:
                return points_to_num[round(point,9)]
            except KeyError:
                return None

        borders = []
        for i in range(init,end):
            trajectory_in_grid = []
            parallel = array(blade.trajectories[i])
            p1, sorted_parallel1 = self._closest_meridian_point(meridian1, parallel, blade)
            p2, sorted_parallel2 = self._closest_meridian_point(meridian2, parallel, blade)
            parallel1 = mathtools.direction_in_halfplane(parallel,p1[3:6])
            parallel2 = mathtools.direction_in_halfplane(parallel,p2[3:6])
            p1, sorted_parallel1 = self._closest_meridian_point(meridian1, parallel1, blade)
            p2, sorted_parallel2 = self._closest_meridian_point(meridian2, parallel2, blade)
            index_left = self._get_first_left_meridian_point_index(parallel, sorted_parallel1, p1)
            index_right = self._get_first_right_meridian_point_index(parallel, sorted_parallel2, p2)

            if index_left==None or index_right==None:
                pass
            elif abs(index_right - index_left)%(len(parallel)-1) == 1:
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

    def _closest_meridian_point(self, meridian, parallel):
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
        blade = self.load_blade()
        model = blade.select_model(closest_meridian_point)
        
        return self._get_ray(model,closest_meridian_point), sorted_parallel

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

    def _get_ray(self, model,point):
        df = model.df(point)
        df = df/linalg.norm(df)
        return array(list(point)+list(df))

    def compute_rays_from_parallels(self, parallels, borders = None):
        """
        This method gets a list of point numbers (parallel db format) and
        create a list of rays (point-normal) with possible border points.
        It removes empty borders and empty parallels, but if both are empty,
        it will return an empty list.

        Keyword arguments:
        parallels -- list of [ list of (point numbers)]
        borders -- If present, must be a list shape (N,2) with begin and end of each parallel
        """

        rays = []
        ntp = self.get_sorted_points()
        blade = self.load_blade()
        parallels = copy.deepcopy(parallels)
        
        db = self.load_db()
        for parallel in parallels:
            for point in parallel:
                if point not in db.keys():
                    parallel.remove(point)
       
        if borders is None:
            for parallel in parallels:
                if len(parallel)==0:
                    rays+=[[]]
                    continue
                model = blade.select_model(ntp[parallel[0]])
                rays += [ map( lambda x: self._get_ray(model,ntp[x]), parallel ) ]

        else:
            for i in range(len(parallels)):
                traj = []
                if len(borders[i][0])>0:
                    model = blade.select_model(borders[i][0])
                    traj.append(self._get_ray(model,borders[i][0]))

                if len(parallels[i])>0:
                    model = blade.select_model(ntp[parallels[i][0]])
                    traj += map( lambda x: self._get_ray(model,ntp[x]), parallels[i])
                    
                if len(borders[i][1])>0:
                    model = blade.select_model(borders[i][1])
                    traj.append(self._get_ray(model,borders[i][1]))
                    
                rays.append(traj)
        return rays
                
    def bases_validation(self, parallels, bases):
        """
        Validate bases for given parallels.
        
        Keyword arguments:
        parallels -- list of rays (N,6)
        bases -- [PSAlpha] (N,3), list of tuples
        """
        
        feasible_bases = []
        blade = self.load_blade()
        for base in bases:
            rp = rail_place.RailPlace(base)
            self.turb.place_rail(rp)
            self.turb.place_robot(rp)
            for parallel in parallels:
                joint_solutions = planning.compute_robot_joints(
                    self.turb, parallel, 0, blade.trajectory_iter_surface)
                if len(joint_solutions) != len(parallel):
                    break
            else:
                feasible_bases.append(base)
            continue
        return feasible_bases

    def remove_point(self, points):
        """
        Remove given points from db, and grid_to_trajectories.
        
        Keyword arguments:
        points (or rays) -- list of points to be removed
        """
        
        db = self.load_db()
        points_to_num = self.load_points_to_num()
        
        for point in points:
            key_point = tuple(point[0:3])     
            try:
                key_num = points_to_num[round(key_point,9)]
            except KeyError:
                continue

            db.pop(key_num,None)
    
        try:
            save_pickle(db, join(self.path,'db.pkl'))
        except IOError: None

        return
                
    def clear_db(self):
        """
        Clear the db.
        """

        db = self.load_db()
        for key, value in db.iteritems():
            db[key] = set()
        try:
            save_pickle(db, join(self.path,'db.pkl'))
        except IOError: raise
        return db

    def clear_visited_bases(self):
        """
        Clear the visited_bases.
        """

        db = self.load_visited_bases()
        for key in db:
            db[key] = False
        try:
            save_pickle(db, join(self.db_main_path,'visited_bases.pkl'))
        except IOError: raise
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

    def plot_bases_db(self, vis):
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
            vis.plot(rp.getXYZ(self.turb.config),'base',(0,0,1))
        return

    def plot_grid(self, grid_num, vis):
        grid_to_trajectories = self.load_grid_to_trajectories()
        blade = self.load_blade()
        trajectories, borders = grid_to_trajectories[grid_num]
        rays = self.compute_rays_from_parallels(blade, trajectories, borders)
        vis.plot_lists(rays, 'rays', color=(1,0,0))
        return
        
