from numpy import *
import coating
from openravepy import *
from openravepy.misc import SpaceSamplerExtra
#====================================================================================================================
env=Environment()
env.Load("../Turbina/env_mh12_0_16.xml")
robot = env.GetRobots()[0]
target = env.GetBodies()[0]
manip = robot.GetActiveManipulator()
handles=[]

ikmodel = databases.inversekinematics.InverseKinematicsModel(robot=robot,iktype=IkParameterization.Type.Transform6D)
if not ikmodel.load():
    ikmodel.autogenerate()

# PARAMETERS
facevector = [1,0,0]
theta = [0,0,0]
coatingdistance = 0.23 # coating distance
coatingdistancetolerance = 0.01
numberofangles = 8 # degree step
tolerance = 30 # degrees
#====================================================================================================================
def solution(pN, ray):
# inputs: - pN (6x1) is the manipulator's base position and orientation.
# - ray (6x1) is the point to be coated
# solution finds joints solutions for a specific point given maximum tolerance.
    _, iksolList, indexlist1 = coating.WorkspaceOnPose(pN, 0, [ray],robot,ikmodel,facevector,theta,coatingdistancetolerance)
    print 'solution: iksolList - ',array(iksolList).shape
    _, AlliksolList, _ = coating.AllExtraCoating2([ray],indexlist1,coatingdistance,numberofangles,tolerance,ikmodel,facevector,coatingdistancetolerance)
    print 'solution: AlliksolList - ',array(AlliksolList).shape
    if iksolList: return iksolList
    elif AlliksolList: return AlliksolList
    else: return False

def fullSolution(pN,ray):
# inputs: - pN (6x1) is the manipulator's base position and orientation.
# - ray (6x1) is the point to be coated
# fullSolution finds all joints solutions for a specific point given tolerance.
    angles = arange(0,tolerance,0.1)
    numberofangles = 10
    for angle in angles:
        angle=1.0*pi*angle/180
        Rv3tol = coating.RotationCalc2([0,0,1],angle)
        p = dot(facevector,transpose(Rv3tol))
        k = 1.0*2*pi/numberofangles
        for i in range(0,numberofangles):
            alfa = k*i
            Rv1alfa = coating.RotationCalc2(facevector,alfa)
            tempP = dot(p,transpose(Rv1alfa))
            _, iksolList, _ = coating.WorkspaceOnPose(pN, 0, [ray],robot,ikmodel,tempP,theta,coatingdistancetolerance)    
            if iksolList:
                return iksolList
    return []

def Qvector_backtrack(y,Q):
# inputs: - y is array(), computed coated points.
# - Q is array(), computed robot's joints.
# Qvector_backtrack method is called when the optimization failed and the previous
# solution should be recomputed.
    for rev_y in reversed(y):
        res = coating.optmizeQ(robot,ikmodel,manip,rev_y,Q[-1])
        Q.append(res.x)
        if not res.success:
            print 'Qvector_backtrack error'
    return list(reversed(Q))

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

def optmizeQ(robot,ikmodel,manip,P,q0):
# inputs: - P 6x1 vector is the goal (point to be coated).
# - q0 is the inicial configuration of the robot (robot's joints).
# optimzeQ minimizes the orientation error.   
    n = [P[3],P[4],P[5]]; P=[P[0],P[1],P[2]]
    def func(q):
        robot.SetDOFValues(q,ikmodel.manip.GetArmIndices())
        T=manip.GetTransform()
        Rx = T[0:3,0]/sqrt(dot(T[0:3,0],T[0:3,0]))
        return -dot(n,Rx)
    def consfunc(q):
        robot.SetDOFValues(q,ikmodel.manip.GetArmIndices())
        T=manip.GetTransform()
        pos = T[0:3,3]
        v = pos-P
        return dot(v,v)
    cons = ({'type':'eq',
             'fun': consfunc})
    res = minimize(func, q0, constraints=cons, method='SLSQP', options={'disp': False})
    return res

def doPath(pN, trajectories):
# inputs: - pN 6x1 vector is the manipulator's base
# position and orientation (x,y,z,nx,ny,nz).
# - trajectories are sorted array(array()) (points to be coated in the right order)
# doPath is the method which iterates points of the trajectories,
# computing optimal robot's joints (minimizing orientation error).
    global handles
    QA = []
    q0=solution(pN,trajectories[0][0])
    q0=q0[0][0]
    robot.SetDOFValues(q0,ikmodel.manip.GetArmIndices())
    for trajectory in trajectories:
        Q=[]
        for index in range(0,len(trajectory)):
            res = optmizeQ(robot,ikmodel,manip,trajectory[index],q0)
            if not coating.CheckDOFLimits(robot,res.x):
                #TODO check collision and angle (tolerance)
                iksolList = fullSolution(pN,trajectory[index])
                Q=[iksolList[0][0]]
                Q = Qvector_backtrack(trajectory[:index],Q)
            else:
                if res.success:
                    Q.append(res.x)
                    q0=res.x
                    handles=coating.plotPoint(env,trajectory[index], handles,array((1,0,0)))
                else: break
        QA.append(Q)
        savez_compressed('trajectory/'+'Q.npz', array=QA)
    return QA
    
def main(pN, BladePosition, trajectories):
    global handles
    coating.RotateBodies(env, BladePosition)
    trajectories = sortTrajectories(pN[0], trajectories)
    savez_compressed('trajectory/'+'croppedTraj.npz', array=trajectories)
    
    Q = doPath(pN, trajectories)
    return Q, trajectories
