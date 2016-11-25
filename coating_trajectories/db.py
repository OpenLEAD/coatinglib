from path_filters import filter_trajectories, side_filter
import planning
from numpy import save, load, random, array
from os.path import basename, splitext, join, exists, isfile
from os import makedirs, listdir
import copy
import rail_place
from colorsys import hls_to_rgb
import cPickle


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

    def __init__(self, path, blade=None, db = None):

        self.path = path
        self.create_db(blade, db)
        del blade
        
        if not exists(self.path):
            makedirs(self.path)
       
    def save_db_pickle(self, obj, name ):
        with open(name, 'wb') as f:
            cPickle.dump(obj, f, cPickle.HIGHEST_PROTOCOL)
        return

    def load_db_pickle(self, name):
        with open(name, 'rb') as f:
            return cPickle.load(f)

    def load_db_npy(self, path):
        return load(path).item()

    def load_db(self):
        """
        Load main database num:num. First num represents the points and
        the second represents the bases.
        """

        name = join(self.path,'fixed_db','db.pkl')
        return self.load_db_pickle(name)

    def create_db(self, blade, db_in):
        """
        Create the database.
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
           
        if db_in is not None:
            db = self.load_db()
            db = self.merge_db(db_in, db)
            self.save_db_pickle(db, join(self.path,'fixed_db','db.pkl'))
            self.save_db_pickle(db_points_to_num, join(self.path,'fixed_db','db_points_to_num.pkl'))
        return

    def load_db_points_to_num(self):
        """
        Load points_to_num database points:num. Points are tuples (x,y,z) and
        the num maps the points (counter to reduce database complexity).

        keyword arguments:
        """

        path = join(self.path, 'fixed_db', 'db_points_to_num.pkl')
        return self.load_db_pickle(path)

    def load_db_bases_to_num(self):
        """
        Load bases_to_num database bases:num. Bases are tuples railplace and
        the num maps the bases (counter to reduce database complexity).

        keyword arguments:
        path -- (optional) The path of the db, if it exists.
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

    def get_bases(self, db):
        bases = set()
        for val in db.values():
            bases = val|bases
        return bases

    def convert_db_point_base(self, db_file, db_points_to_num = None):
        """
        Method to convert a point:base dict file to dict num:num.

        Keyword arguments:
        db_file -- dict point:base
        """

        converted_db = dict()

        if db_points_to_num is None:
            db_points_to_num = self.load_db_points_to_num()

        db_bases_to_num = self.load_db_bases_to_num() 
        bases = self.get_bases(db_file)
        for i in range(0,len(bases)):
            base = bases.pop()
            db_bases_to_num[base] = db_bases_to_num.get(base,set([len(db_bases_to_num)]))
        self.save_db_pickle(db_bases_to_num, join(self.path,'fixed_db','db_bases_to_num.pkl'))
        
        for key, values in db_file.iteritems():
            intkey = db_points_to_num[key]
            for value in values:
                converted_db[intkey] = db_bases_to_num[value]
        return converted_db

    def convert_db_point_base_directory(self, path, directory_to_save):
        """
        Method to convert all databases in a specific folder and
        save the converted dbs to a given folder.

        keyword arguments:
        path -- where the databases are
        directory_to_save -- where to save the converted databases
        """

        if not exists(directory_to_save):
            makedirs(directory_to_save)

        db_points_to_num = self.load_db_points_to_num()
        onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
         
        for afile in onlyfiles:
            filename, file_extension = splitext(afile)
            if file_extension == '.npy':
                db_file = self.load_db_npy(join(path,afile))
            elif file_extension == '.pkl':
                db_file = self.load_db_pickle(join(path,afile))
            else:
                continue
            print "Converting file: ", afile
            self.save_db_pickle(self.convert_db_point_base(db_file, db_points_to_num),
                                join(directory_to_save, filename+'.pkl'))
        return


    def merge_db(self, db_file, main_db):
        """
        Merge new db_file to the main database.

        keyword arguments:
        db_file -- file to be merged;
        db (optional) -- main database file;
        """

        for key, value in db_file.iteritems():
            main_db[key] = main_db.get(key,set()) | db_file.get(key)
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
         
        for afile in onlyfiles:
            filename, file_extension = splitext(afile)
            if file_extension == '.pkl':
                print "merging file: ", filename
                db = self.merge_db(self.load_db_pickle(join(path,afile)), db)
        self.save_db_pickle(db, join(self.path,'fixed_db','db.pkl'))
        return

    def generate_db(self, turbine, blade, rail_position, minimal_number_of_points_per_trajectory = 2): 
        """
        Function that store feasible trajectories for random rail_positions. Each feasible point
        in a trajectory is mapped to a rail position that can coat it. Therefore, the DB returns
        a dictionary where points are keys and rail positions are values, e.g.:
             {     (1.0, 1.1, 1.3):                   (2.4   ,  1.2  ), ( 2.4  ,  1.3  )}
             {point( x ,  y ,  z ): rail poisitions   (rail_1, rail_2), (rail_1, rail_2)}

        Keyword arguments:
        turbine -- turbine object
        rail_position -- number of random rail positions
        minimal_number_of_points_per_trajectory -- trajectories smaller then
        this number of points are croped
        """

        db_points_to_num = self.load_db_points_to_num()
        db_bases_to_num = self.load_db_bases_to_num()
        
        db = dict()
        base = tuple(rail_position.getPSAlpha())
        db_bases_to_num[base] = db_bases_to_num.get(base, set([len(db_bases_to_num)]))
        self.save_db_pickle(db_points_to_num, join(self.path,'fixed_db','db_bases_to_num.pkl'))

        turbine.place_rail(rail_position)
        turbine.place_robot(rail_position)
            
        if turbine.check_rail_collision():
            return
        if turbine.check_robotbase_collision():
            return

        filtered_trajectories = filter_trajectories(turbine, blade.trajectories)

        for filtered_trajectory in filtered_trajectories:
            for filtered_trajectory_part in filtered_trajectory:
                evaluated_points = 0
                while evaluated_points < len(filtered_trajectory_part):
                    try:
                        lower, _, _ = planning.compute_first_feasible_point(turbine, filtered_trajectory_part[evaluated_points:], blade.trajectory_iter_surface)
                        evaluated_points = evaluated_points + lower
                    except ValueError: 
                        evaluated_points = len(filtered_trajectory_part)
                        continue

                    joint_solutions = planning.compute_robot_joints(turbine, filtered_trajectory_part, evaluated_points, blade.trajectory_iter_surface)
                    upper = evaluated_points + len(joint_solutions)

                    if upper-evaluated_points > minimal_number_of_points_per_trajectory:
                        for point in filtered_trajectory_part[evaluated_points:upper]:
                            db[db_points_to_num[tuple(point[0:3])]] = db_bases_to_num[base]
                        print 'saved'
                    evaluated_points = upper+1
        return db

    def plot_points_gradient(self, vis):
        N = 0
        db = self.load_db()
        for key, value in db.iteritems():
            N = max(N,len(value))

        db_points_to_num = self.load_db_points_to_num()
        all_points_tuple, all_points_num = db_points_to_num.keys(), db_points_to_num.values()
        del db_points_to_num
        all_points_tuple = [x for (y,x) in sorted(zip(all_points_num,all_points_tuple))]
        del all_points_num
        
        index = map(int,random.uniform(0,len(db)-1,int(len(db)*1.0/200)))
        points_num = array(db.keys())[index]
        bases_num = array(db.values())[index]

        for i in range(0,len(points_num)):
            vis.plot(all_points_tuple[points_num[i]], 'points_gradient',
                     color = hls_to_rgb(len(bases_num[i])*1.0/(3*N),0.5,1))
        return

    def plot_points_db(self, db, vis):
        index = map(int,random.uniform(0,len(db)-1,int(len(db)*1.0/100)))
        points_num = array(db.keys())[index]

        db_points_to_num = self.load_db_points_to_num()
        all_points_tuple, all_points_num = db_points_to_num.keys(), db_points_to_num.values()
        del db_points_to_num
        all_points_tuple = [x for (y,x) in sorted(zip(all_points_num,all_points_tuple))]
        del all_points_num
        
        for i in points_num:
            vis.plot(all_points_tuple[i], 'points_db')
        return

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
