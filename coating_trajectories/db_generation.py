from path_filters import filter_trajectories
import planning
from numpy import save
from visualizer import Visualizer


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
    vis = Visualizer(turbine.env)
    for rp in rail_positions:
        turbine.place_rail(rp)
        turbine.place_robot(rp)
        
        if turbine.check_rail_collision():
            continue
        
        if turbine.check_robotbase_collision():
            continue

        filtered_trajectories = filter_trajectories(turbine, blade.trajectories)

##        vis.remove_points('trajectories')
##        for filtered_trajectory in filtered_trajectories:
##            for filtered_trajectory_part in filtered_trajectory:
##                vis.plot(filtered_trajectory_part, 'trajectories', ((0,0,1)))
##        x = raw_input('wait')

        counter = 0
        for filtered_trajectory in filtered_trajectories:
            for filtered_trajectory_part in filtered_trajectory:
                try:
                    lower, _, _ = planning.compute_first_feasible_point(turbine, filtered_trajectory_part)
                except ValueError: #raise
                    continue

                joint_solutions = planning.compute_robot_joints(turbine, filtered_trajectory_part, lower)
                upper = lower + len(joint_solutions)
                
##                for i in range(0,len(joint_solutions)):
##                    w, alpha = planning.compute_angular_velocities(turbine, joint_solutions, i)
##                    torques = planning.torque_computation(turbine, joint_solutions[i], w, alpha)
##                    velocity_tan_error, position_normal_error,
##                    position_perp_error, angle_error = planning.sensibility(turbine, filtered_trajectory_part[lower+i], w, alpha)

                if upper-lower>minimal_number_of_points_per_trajectory:
                    for point in filtered_trajectory_part[lower:upper]:
                        DB_dict[tuple(point[0:3])] = (DB_dict.get(tuple(point[0:3]),set()) |
                                                      set([tuple(rp.getPSAlpha())]))

            print counter
            save('db.npy', DB_dict) 
            counter += 1
                                                          
    return DB_dict

def save_obj(obj, name ):
    with open('obj/'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open('obj/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)

        
