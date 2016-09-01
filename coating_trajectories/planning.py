from numpy import sqrt, dot, concatenate, arange, array
from openravepy import IkFilterOptions
from math import pi, cos, sin, atan2
from scipy.optimize import minimize

""" Main package for robot joints' positions and velocities planning,
robot base calculation, torque and manipulability analysis.

Keyword arguments:
places - robot places to coat a specific set of parallels
turbine - full environment and parameters
"""
_filter_options = {'mh12': mh12_filter}

def _std_robot_filter(turbine, trajectories):
    raise ValueError("No trajectory filter for "+turbine.robot.GetName()+" robot. Create new function.")

def mh12_filter(turbine, trajectories):
    def distance_robot(point,turbine):
        pass # TODO
    filtered_trajectories = []
    for jectory in trajectories:
        
        filtered_jectory = []
        for point in jectory:
            
            if distance_robot(point[0:3],turbine):
                pass #TODO
    

def filter_trajectories(turbine, trajectories):
    name = turbine.robot.GetName()
    _filter_options.get(name,_std_robot_filter)(turbine, trajectories)

def sort_trajectories(trajectories):
    """ Arrange the trajectories in a zigzagging way.
    """
    sorted_trajectories = []
    for trajectory in trajectories:
        sorted_trajectory = []

def sortTrajectories2(trajectories):
    """ Arrange the trajectories in a zigzagging way.
    """
    sortedTrajectories = []
    i=1
    for trajectory in trajectories:
        if len(trajectory)>1:
            theta = []
            for point in trajectory:
                theta.append(dot(point[0:3], [0,1,0]))
            theta = array(theta)
            sortedTrajectory = []
            sortedTrajectory.extend(trajectory[theta.argsort()])
           # if i%2:
           #     sortedTrajectory.reverse()
            sortedTrajectories.append(sortedTrajectory)
        elif len(trajectory)==1:
            sortedTrajectories.append(trajectory)
        i+=1    
    return array(sortedTrajectories)

  
def sortTrajectories(x, trajectories):
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
                theta.append(atan2(-point[2],point[0]))
            theta = array(theta)
            sortedTrajectory = []
            if len(theta[theta>0])>0:    
                sortedTrajectory.extend((trajectory[theta>0])[(theta[theta>0]).argsort()])
            if len(theta[theta<0])>0:
                if x<0:
                    sortedTrajectory.extend((trajectory[theta<0])[(theta[theta<0]).argsort()])
                else:
                    En = (trajectory[theta<0])[(theta[theta<0]).argsort()]
                    En.extend(sortedTrajectory)
                    sortedTrajectory = En
            if i%2:
                sortedTrajectory.reverse()
            sortedTrajectories.append(sortedTrajectory)
        elif len(trajectory)==1:
            sortedTrajectories.append(trajectory)
        i+=1    
    return sortedTrajectories

def doPath(place, turbine):
    """ Iterates points of the trajectories, computing optimal robot's joints
    (minimizing orientation error).

    Keyword arguments:
    x - x position of the robot
    """

    turbine.robot.GetLink('Flame').Enable(False)
    QA = []
    q0, _ = InverseKinematic(place.T, turbine, place.trajectories[0][0])
    q0 = q0[0][0]
    turbine.robot.SetDOFValues(q0, turbine.ikmodel.manip.GetArmIndices())

    for trajectory in place.trajectories:
        Q=[]
        for index in range(0,len(trajectory)):
            res = optmizeQ(trajectory[index], q0)
            if not CheckDOFLimits(turbine.robot, res.x):
                #TODO check collision and angle (tolerance)
                iksolList = InverseKinematic(place.T, turbine, trajectory[index])
                Q=[iksolList[0][0]]
                Q = Qvector_backtrack(trajectory[:index],Q)
            else:
                if res.success:
                    Q.append(res.x)
                    q0 = res.x
                else: break
        QA.append(Q)
    return QA

def compute_joints(places):
    for place in places:
        x = place.T[3,0] 
        place.trajectories = sortTrajectories(x, place.trajectories)
        Q = doPath(place, turbine)

def optmizeQ(turbine, P, q0):
    """ Minimizes the orientation error.

    Keyword arguments:
    P -- 6x1 vector is the goal (point to be coated).
    q0 -- is the inicial configuration of the robot (robot's joints).
    """
    n = [P[3],P[4],P[5]]; P=[P[0],P[1],P[2]]
    def func(q):
        turbine.robot.SetDOFValues(q, turbine.ikmodel.manip.GetArmIndices())
        T = turbine.robot.GetActiveManipulator().GetTransform()
        Rx = T[0:3,0]/sqrt(dot(T[0:3,0],T[0:3,0]))
        return -dot(n,Rx)
    def consfunc(q):
        turbine.robot.SetDOFValues(q, turbine.ikmodel.manip.GetArmIndices())
        T = turbine.robot.GetActiveManipulator().GetTransform()
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
        res = coating.optmizeQ(turbine, rev_y, Q[-1])
        Q.append(res.x)
        if not res.success:
            print 'Qvector_backtrack error'
    return list(reversed(Q))

def InverseKinematic(T, turbine, point):
    """ Pose the robot and solve inverse kinematics given point (IKFast).
    """
    
    facevector = [1,0,0] # robot end-effector direction 
    turbine.robot.SetTransform(T)
    coating_tolerance = turbine.coating.max_distance-turbine.coating.ideal_distance

    iksol = IKFast(turbine.ikmodel, facevector, point)
    if len(iksol)>0:
        return iksol, True
 
    if coating_tolerance!=0:
        iksol = IKFast(turbine.ikmodel,
                       facevector,
                       concatenate((point[0:3]+coating_tolerance*point[3:6],
                                    point[3:6]))
                       )
        if len(iksol)>0:
            return iksol, True
        
    if turbune.coating.angle_tolerance>0:
        angles = arange(0, turbune.coating.angle_tolerance, 0.001)
        numberofangles = 10
        for angle in angles:
            angle=1.0*pi*angle/180
            Rv3tol = Raxis([0,0,1],angle)
            p = dot(facevector,transpose(Rv3tol))
            k = 2*pi/numberofangles
            for i in range(0,numberofangles):
                alfa = k*i
                Rv1alfa = coating.Raxis(facevector,alfa)
                iksol = IKFast(turbine.ikmodel,
                               facevector,
                               dot(p,
                                   transpose(Rv1alfa))
                               )
                if len(iksol)>0:
                    return iksol, True
            else: continue
            break     
        else: 
            return [], False
            
def IKFast(ikmodel, facevector, point):
    """ Call openrave IKFast.
    """

    T = zeros((4,4))
    T[0:3,0:3] = Rab(facevector, point[3:6])
    T[0:3,3] = point[0:3]
    T[3,0:4] = [0,0,0,1]
    iksol = ikmodel.manip.FindIKSolutions(T, IkFilterOptions.CheckEnvCollisions)
        
    return iksol

def CheckDOFLimits(robot, Q):
    Joints = robot.GetJoints()
    for i in range(0,len(Q)):
        l,u = Joints[i].GetLimits()
        l = l[0]
        u = u[0]
        if not Q[i]>=l+0.1 or not Q[i]<=u-0.1:
            return False
    return True    
