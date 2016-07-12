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

manipulabilityNUM = 0
globalManipPos = []
globalManipOri = []

def solution(pN, ray):
    _, iksolList, indexlist1 = coating.WorkspaceOnPose(pN, 0, [ray],robot,ikmodel,facevector,theta,coatingdistancetolerance)
    print 'solution: iksolList - ',array(iksolList).shape
    _, AlliksolList, _ = coating.AllExtraCoating2([ray],indexlist1,coatingdistance,numberofangles,tolerance,ikmodel,facevector,coatingdistancetolerance)
    print 'solution: AlliksolList - ',array(AlliksolList).shape
    if iksolList: return iksolList
    elif AlliksolList: return AlliksolList
    else: return False

def fullSolution(ray):
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
    for rev_y in reversed(y):
        res = coating.optmizeQ(robot,ikmodel,manip,rev_y,Q[-1])
        Q.append(res.x)
        if not res.success:
            print 'Qvector_backtrack error'
    return list(reversed(Q))

def isViable(q0,norm):
    # Alterar
    global manipulabilityNUM
    global globalManipPos
    global globalManipOri
    manipulability, manipjpos, manipjori = coating.manipulabilityS(manip)
    globalManipPos.append(manipjpos)
    globalManipOri.append(manipjori)

    manipjpos = 0.3*manipjpos + 0.7*manipulabilityNUM
    if manipulabilityNUM > manipjpos and manipulabilityNUM<=0.15:
        print "LOW MANIPULABILITY: ", manipjpos
        return False
    manipulabilityNUM = manipjpos
    
    if len(q0)>0:
        robot.SetDOFValues(q0,ikmodel.manip.GetArmIndices())
        T=manip.GetTransform()
        Rx = T[0:3,0]/sqrt(dot(T[0:3,0],T[0:3,0]))
        if dot(norm,Rx) >= math.cos(tolerance*math.pi/180):
            return True
        else:
            print "ANGLE TOLERANCE FAIL: ", 180*math.acos(dot(norm,Rx))/math.pi
            return False
    else:
        return False

def sortTrajectories(pNx, trajectories):
    sortedTrajectories = []
    i=0
    for trajectory in trajectories:
        if len(trajectory)>1:
            theta = []
            for point in trajectory:
                theta.append(math.atan2(-point[2],point[0]))
            theta=array(theta)
            St = []; En = []
            if len(theta[theta>0])>0:    
                St = [x for (y,x) in sorted(zip(theta[theta>0],trajectory[theta>0]))]
            if len(theta[theta<0])>0:
                En = [x for (y,x) in sorted(zip(theta[theta<0],trajectory[theta<0]))]
            if St and En:
                if pNx<0:
                    print 'St = ', St
                    print 'En = ', En
                    sortedTrajectory=concatenate((St,En))
                else:    
                    sortedTrajectory=concatenate((St,En))
            if i%2:
                sortedTrajectory.reverse()
            sortedTrajectories.append(sortedTrajectory)
        elif len(trajectory)==1:
            sortedTrajectories.append(trajectory)
    return sortedTrajectories        

def doPath(pN, trajectories):
    QA = []
    q0=solution(trajectories[0][0])
    q0=q0[0][0]
    robot.SetDOFValues(q0,ikmodel.manip.GetArmIndices())
    for trajectory in trajectories:
        Q=[q0]
        for index in range(0,len(trajectory[1:])):
            res = coating.optmizeQ(robot,ikmodel,manip,trajectory[index],q0)
            if not coating.CheckDOFLimits(robot,res.x):
                #TODO Verificar colisao tb
                iksolList = fullSolution(point)
                Q=[iksolList[0][0]]
                Q = Qvector_backtrack(trajectory[:index],Q)
            else:
                if res.success:
                    Q.append(res.x)
                    q0=res.x
                    handles=coating.plotPoint(point, handles,array((0,1,0)))
                else: break
        QA.append(Q)        
    return QA
    
def main(pN, BladePosition, trajectories):
    global handles
    coating.RotateBodies(env, BladePosition)
    trajectories = sortTrajectories(pN[0], trajectories)
    
    env.SetViewer('qtcoin')
    Q = doPath(pN, trajectories)
    return Q

def realCoatedPoints(Q):
    global handles
    newY=[]
    for q in Q:
        robot.SetDOFValues(q,ikmodel.manip.GetArmIndices())
        T=manip.GetTransform()
        Rx = array(T[0:3,0]/sqrt(dot(T[0:3,0],T[0:3,0])))
        p = T[0:3,3]
        newY.append(p-0.23*Rx)
    handles=plotPoints(newY, handles,array((0,0,1)))
    return handles

def getPointsfromQ(Q):
    points=[]
    for q in Q:
        robot.SetDOFValues(q,ikmodel.manip.GetArmIndices())
        T=manip.GetTransform()
        points.append(T[0:3,3])
    return array(points)  
