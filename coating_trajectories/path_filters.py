from numpy import array, sum, nonzero, concatenate, split, dot
from mathtools import direction_in_halfplane


def _std_robot_filter(turbine, *args):
    raise ValueError("No trajectory filter for " + turbine.robot.GetName() + " robot. Create new function.")

def _mh12_points_filter(turbine, points):
    pistol = turbine.robot.GetLinks()[6].GetGeometries()[0].GetCylinderHeight()
    flame = turbine.robot.GetLinks()[7].GetGeometries()[0].GetCylinderHeight()
    _working_radius_squared = (1.285 + pistol + flame) ** 2  ##HARDCODED
    _non_working_radius_squared = (0.410) ** 2

    def distance_robot_squared(trajectory_points):
        delta = turbine.robot.GetJoints()[1].GetAnchor() - trajectory_points
        return sum(delta * delta, 1)

    points = array(points)
    distances = distance_robot_squared(points[:, 0:3])
    boo = (distances < _working_radius_squared) & (distances > _non_working_radius_squared)
    return points[boo]


def points_side_filter(turbine, points):
    Rx = -turbine.robot.GetTransform()[0:3, 0]
    return direction_in_halfplane(points, Rx)

_points_filter_options = {'mh12': _mh12_points_filter}


def filter_points(turbine, points, do_side_filter=True):
    name = turbine.robot.GetName()
    if do_side_filter:
        return _points_filter_options.get(name, _std_robot_filter)(turbine, points_side_filter(turbine, points))
    else:
        return _points_filter_options.get(name, _std_robot_filter)(turbine, points)


##### TRAJECTORY FILTERS #####


def _mh12_trajectory_filter(turbine, trajectories, N):
    pistol = turbine.robot.GetLinks()[6].GetGeometries()[0].GetCylinderHeight()
    flame = turbine.robot.GetLinks()[7].GetGeometries()[0].GetCylinderHeight()
    _working_radius_squared = (1.285 + pistol + flame) ** 2  ##HARDCODED
    _non_working_radius_squared = (0.410) ** 2

    def distance_robot_squared(trajectory_points):
        delta = turbine.robot.GetJoints()[1].GetAnchor() - trajectory_points
        return sum(delta * delta, 1)

    filtered_trajectories = []
    for trajectory in trajectories:
        trajectory = array(trajectory)
        distances = distance_robot_squared(trajectory[:, 0:3])
        boo = (distances < _working_radius_squared) & (distances > _non_working_radius_squared)
        indices = nonzero(boo[1:] != boo[:-1])[0] + 1
        split_trajectory = split(trajectory, indices)
        filtered_trajectory_list = split_trajectory[0::2] if boo[0] else split_trajectory[1::2]
        if boo[0] & boo[-1] & (len(filtered_trajectory_list) > 1):
            filtered_trajectory_list[0] = concatenate((filtered_trajectory_list[-1], filtered_trajectory_list[0]))
            del filtered_trajectory_list[-1]

        filtered_trajectory_list = filter(lambda x: len(x) > N, filtered_trajectory_list)
        if len(filtered_trajectory_list) > 0:
            filtered_trajectories.append(filtered_trajectory_list)

    return filtered_trajectories


def trajectory_side_filter(turbine, trajectories):
    Rx = -turbine.robot.GetTransform()[0:3, 0]

    filtered_trajectories = []
    for trajectory in trajectories:
        filtered_trajectory = direction_in_halfplane(trajectory, Rx)
        if len(filtered_trajectory) > 0:
            filtered_trajectories += [filtered_trajectory]

    return filtered_trajectories


_trajectory_filter_options = {'mh12': _mh12_trajectory_filter}


def filter_trajectories(turbine, trajectories, N=100, do_side_filter=True):
    name = turbine.robot.GetName()
    if do_side_filter:
        return _trajectory_filter_options.get(name, _std_robot_filter)(turbine,
                                                                       trajectory_side_filter(turbine, trajectories), N)
    else:
        return _trajectory_filter_options.get(name, _std_robot_filter)(turbine, trajectories, N)


