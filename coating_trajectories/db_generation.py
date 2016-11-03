from path_filters import filter_trajectories
import planning
from numpy import save, load
import logging
import time
from os.path import basename, splitext, join, exists, isfile
from datetime import datetime
from os import makedirs, listdir
import copy

def generate_db(turbine, blade, rail_positions, DB_dict = dict(), minimal_number_of_points_per_trajectory = 100): 
    """
    Function that store feasible trajectories for random rail_positions. Each feasible point
    in a trajectory is mapped to a rail position that can coat it. Therefore, the DB returns
    a dictionary where points are keys and rail positions are values, e.g.:
         {     (1.0, 1.1, 1.3):                   (2.4   ,  1.2  ), ( 2.4  ,  1.3  )}
         {point( x ,  y ,  z ): rail poisitions   (rail_1, rail_2), (rail_1, rail_2)}

    Keyword arguments:
    turb -- turbine object
    blade -- blade object
    rail_positions -- number of random rail positions
    """
    
    now = datetime.now()
    _ , week, dayofweek = now.isocalendar()
    directory = join('week'+str(week), 'day'+str(dayofweek))
    if not exists(directory):
        makedirs(directory)

    logging.basicConfig(filename=join(directory,'trajectory_constraints' + now.strftime('%X').replace(':','_') + '.log'), level=logging.DEBUG)
    for rp in rail_positions:
        turbine.place_rail(rp)
        turbine.place_robot(rp)
        
        if turbine.check_rail_collision():
            logging.info('Bad rail placement: rail collision')
            continue
        
        if turbine.check_robotbase_collision():
            logging.info('Bad robot placement: robot base collision')
            continue

        filtered_trajectories = filter_trajectories(turbine, blade.trajectories)

        counter = 0
        for filtered_trajectory in filtered_trajectories:
            logging.info('Trajectory: '+str(counter))
            for filtered_trajectory_part in filtered_trajectory:
                evaluated_points = 0
                while evaluated_points < len(filtered_trajectory_part):
                    try:
                        lower, _, _ = planning.compute_first_feasible_point(turbine, filtered_trajectory_part[evaluated_points:], blade.trajectory_iter_surface)
                        evaluated_points = evaluated_points + lower
                    except ValueError: #raise
                        evaluated_points = len(filtered_trajectory_part)
                        continue

                    joint_solutions = planning.compute_robot_joints(turbine, filtered_trajectory_part, evaluated_points, blade.trajectory_iter_surface)
                    upper = evaluated_points + len(joint_solutions)
                    
    ##                for i in range(0,len(joint_solutions)):
    ##                    w, alpha = planning.compute_angular_velocities(turbine, joint_solutions, i)
    ##                    torques = planning.torque_computation(turbine, joint_solutions[i], w, alpha)
    ##                    velocity_tan_error, position_normal_error,
    ##                    position_perp_error, angle_error = planning.sensibility(turbine, filtered_trajectory_part[lower+i], w, alpha)

                    if upper-evaluated_points > minimal_number_of_points_per_trajectory:
                        logging.info('Number of points saved in database: '+str(upper-evaluated_points))
                        for point in filtered_trajectory_part[evaluated_points:upper]:
                            DB_dict[tuple(point[0:3])] = (DB_dict.get(tuple(point[0:3]),set()) |
                                                          set([tuple(rp.getPSAlpha())]))
                    else: logging.info('Points not saved in database, because number of coated points is: '+str(upper-evaluated_points))
                    logging.info('\n')
                    evaluated_points = upper+1
            counter+=1
    return DB_dict

def invert_db(db, parallels, regions):
    # Regions as list of [list of tuple (parallel_index, begin_index, end_index)]
    psa_db = list()
    
    for region in regions:
        inverse_set = set()
        
        for parallel_index, begin_index, end_index in region:
            for point in parallels[parallel_index, begin_index:end_index]:
                inverse_set |= db[point]

        psa_db.append(inverse_set)
            

def save_obj(obj, name ):
    with open('obj/'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open('obj/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)

def save_db(db_file, directory_to_save):
    now = datetime.now()
    _ , week, dayofweek = now.isocalendar()
    directory = join(directory_to_save,'week'+str(week), 'day'+str(dayofweek))
    if not exists(directory):
        makedirs(directory)
    filename=join(directory,'db' + now.strftime('%X').replace(':','_') + '.npy')
    save(filename, db_file)
    return

def load_db(directory_to_load):
    return load(directory_to_load).item()
        
