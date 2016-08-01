from mathtools import InverseKinematic, CheckDOFLimits
from numpy import minimize, sqrt, dot

class Planning:
    """ Main class for robot joints' positions and velocities planning,
    robot base calculation, torque and manipulability analysis.

    Keyword arguments:
    places - robot places to coat a specific set of parallels
    trajectories - set of parallels (array of arrays)
    """

    def __init__(self, places):
        self.places = places
        self.turbine = places.turbine
        
    def sortTrajectories(self, x):
    """ Arrange the trajectories in a zigzagging way.

    Keyword arguments:
    x - x position of the robot
    """
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
                    if x<0:
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

    def doPath(self, place):
    """ Iterates points of the trajectories, computing optimal robot's joints
    (minimizing orientation error).

    Keyword arguments:
    x - x position of the robot
    """
    QA = []
    q0, _ = InverseKinematic(place.T, self.turbine, place.trajectories[0][0])
    q0 = q0[0][0]
    self.turbine.robot.SetDOFValues(q0, self.turbine.robot.ikmodel.manip.GetArmIndices())
    
    for trajectory in trajectories:
        Q=[]
        for index in range(0,len(trajectory)):
            res = optmizeQ(trajectory[index], q0)
            if not CheckDOFLimits(self.turbine.robot, res.x):
                #TODO check collision and angle (tolerance)
                iksolList = InverseKinematic(place.T, self.turbine, trajectory[index])
                Q=[iksolList[0][0]]
                Q = Qvector_backtrack(trajectory[:index],Q)
            else:
                if res.success:
                    Q.append(res.x)
                    q0 = res.x
                else: break
        QA.append(Q)
    return QA

    def compute_joints(self):
        for place in self.places:
            x = place.T[3,0] 
            trajectories = sortTrajectories(x, place.trajectories)
            Q = doPath(place, trajectories)

    def optmizeQ(self, P, q0):
    """ Minimizes the orientation error.

    Keyword arguments:
    P -- 6x1 vector is the goal (point to be coated).
    q0 -- is the inicial configuration of the robot (robot's joints).
    """
        n = [P[3],P[4],P[5]]; P=[P[0],P[1],P[2]]
        def func(q):
            self.turbine.robot.SetDOFValues(q, self.turbine.ikmodel.manip.GetArmIndices())
            T = self.turbine.robot.GetActiveManipulator().GetTransform()
            Rx = T[0:3,0]/sqrt(dot(T[0:3,0],T[0:3,0]))
            return -dot(n,Rx)
        def consfunc(q):
            self.turbine.robot.SetDOFValues(q, self.turbine.ikmodel.manip.GetArmIndices())
            T = self.turbine.robot.GetActiveManipulator().GetTransform()
            pos = T[0:3,3]
            v = pos-P
            return dot(v,v)
        cons = ({'type':'eq',
                 'fun': consfunc})
        res = minimize(func, q0, constraints=cons, method='SLSQP', options={'disp': False})
        return res    

    def Qvector_backtrack(y, Q):
    """ Call when optimization fails. Previous solution should be recomputed.

    Keyword arguments:
    y -- computed coated points.
    Q -- computed robot's joints.
    """
        for rev_y in reversed(y):
            res = coating.optmizeQ(robot,ikmodel,manip,rev_y,Q[-1])
            Q.append(res.x)
            if not res.success:
                print 'Qvector_backtrack error'
        return list(reversed(Q))
