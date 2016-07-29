

class Planning:
    """ Main class for robot joints' positions and velocities planning,
    robot base calculation, torque and manipulability analysis.

    Keyword arguments:
    places - robot places to coat a specific set of parallels
    trajectories - set of parallels (array of arrays)
    """

    def __init__(self, places, trajectories):
        self.places = places
        self._trajectories = trajectories
        self._turbine = self.places._turbine
        self._robot = self.places._robot
        
    def sortTrajectories(pNx, trajectories):
        
    # inputs: - pNx is the manipulator's base x position
    # - trajectories are non-sorted array(array()), points to be coated.
    # sortTrajectories is the method which iterates points of the trajectories,
    # sorting those points in the right order to be coated, as a cropped path
    # should not keep the right order.
        sortedTrajectories = []
        i=1
        for trajectory in trajectories:
            if len(trajectory)>1:
                theta = []
                for point in trajectory:
                    theta.append(math.atan2(-point[2],point[0]))
                theta=array(theta)
                sortedTrajectory = []
                if len(theta[theta>0])>0:    
                    sortedTrajectory.extend([x for (y,x) in sorted(zip(theta[theta>0],trajectory[theta>0]))])
                if len(theta[theta<0])>0:
                    if pNx<0:
                        sortedTrajectory.extend([x for (y,x) in sorted(zip(theta[theta<0],trajectory[theta<0]))])
                    else:
                        En = [x for (y,x) in sorted(zip(theta[theta<0],trajectory[theta<0]))]
                        En.extend(sortedTrajectory)
                        sortedTrajectory = En
                if i%2:
                    sortedTrajectory.reverse()
                sortedTrajectories.append(sortedTrajectory)
            elif len(trajectory)==1:
                sortedTrajectories.append(trajectory)
            i+=1    
        return sortedTrajectories   
