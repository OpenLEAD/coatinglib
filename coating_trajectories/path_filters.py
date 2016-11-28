from numpy import array, sum, nonzero, concatenate, split, dot

def _std_robot_filter(turbine, trajectories):
    raise ValueError("No trajectory filter for "+turbine.robot.GetName()+" robot. Create new function.")

def _mh12_filter(turbine, trajectories, N):
    pistol = 0.3
    flame = 0.23
    _working_radius_squared = (1.285+pistol+flame)**2
    _non_working_radius_squared = (0.410)**2
    
    def distance_robot_squared(trajectory_points):
        delta = turbine.robot.GetJoints()[1].GetAnchor() - trajectory_points
        return sum(delta*delta,1)

    filtered_trajectories = []
    for trajectory in trajectories:
        trajectory = array(trajectory)
        distances = distance_robot_squared(trajectory[:,0:3])
        boo = (distances < _working_radius_squared) & (distances > _non_working_radius_squared)
        indices = nonzero(boo[1:] != boo[:-1])[0] + 1
        split_trajectory = split(trajectory, indices)
        filtered_trajectory_list = split_trajectory[0::2] if boo[0] else split_trajectory[1::2]
        if boo[0] & boo[-1] & (len(filtered_trajectory_list)>1):
            filtered_trajectory_list[0] = concatenate((filtered_trajectory_list[-1],filtered_trajectory_list[0]))
            del filtered_trajectory_list[-1]

        filtered_trajectory_list = filter(lambda x: len(x)>N,filtered_trajectory_list)
        if len(filtered_trajectory_list)>0:
            filtered_trajectories.append(filtered_trajectory_list)

    return filtered_trajectories

def side_filter(turbine, trajectories):
    Rx = turbine.robot.GetTransform()[0:3,0]
    
    filtered_trajectories = []
    for trajectory in trajectories:
        trajectory = array(trajectory)
        filtered_trajectory = trajectory[ (dot(trajectory[:,3:6],Rx)) < 0 ]
        if len(filtered_trajectory) > 0:
            filtered_trajectories.append(filtered_trajectory)

    return filtered_trajectories
            
    

_filter_options = {'mh12': _mh12_filter}    

def filter_trajectories(turbine, trajectories, N = 100):
    name = turbine.robot.GetName()
    return _filter_options.get(name,_std_robot_filter)(turbine, side_filter(turbine,trajectories), N)
