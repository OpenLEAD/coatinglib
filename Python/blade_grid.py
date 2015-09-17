from numpy import *
from openravepy import *
from openravepy.misc import SpaceSamplerExtra
env=Environment()
env.SetViewer('qtcoin')
env.Load("/home/renan/Documents/EMMA/Turbina/env_pa.xml")
target = env.GetBodies()[0]

# PARAMETERS
delta = 0.02
normalanglerange = 0
directiondelta = 0.4
anglerange = pi/6
coatingdistance = 0.23
robottobladedistance = 0.72
initialdistance = 0.23

def PointsToReach(env, target, delta, normalanglerange,directiondelta):
    cc = RaveCreateCollisionChecker(env,'ode')
    if cc is not None:
            ccold = env.GetCollisionChecker()
            env.SetCollisionChecker(cc)
            cc = ccold

    target.SetTransform(eye(4))
    ab = target.ComputeAABB()
    p = ab.pos()
    e = ab.extents()+0.02 # increase since origin of ray should be outside of object
    sides = array(((0,0,e[2],0,0,-1,e[0],0,0,0,e[1],0),
                 (0,0,-e[2],0,0,1,e[0],0,0,0,e[1],0),
                 (0,e[1],0,0,-1,0,e[0],0,0,0,0,e[2]),
                 (0,-e[1],0,0,1,0,e[0],0,0,0,0,e[2]),
                 (e[0],0,0,-1,0,0,0,e[1],0,0,0,e[2]),
                 (-e[0],0,0,1,0,0,0,e[1],0,0,0,e[2])))
    maxlen = 2*sqrt(sum(e**2))+0.03
    approachrays = zeros((0,6))
    side=sides[1];
    ex = sqrt(sum(side[6:9]**2))
    ey = sqrt(sum(side[9:12]**2))
    if ex/delta > 1000:
         raise ValueError('object is way too big for its discretization! %f > 1000'%(ex/delta))
    XX,YY = meshgrid(r_[arange(-ex,-0.25*delta,delta),0,arange(delta,ex,delta)],
                          r_[arange(-ey,-0.25*delta,delta),0,arange(delta,ey,delta)])
    localpos = outer(XX.flatten(),side[6:9]/ex)+outer(YY.flatten(),side[9:12]/ey)
    N = localpos.shape[0]
    rays = c_[tile(p+side[0:3],(N,1))+localpos,maxlen*tile(side[3:6],(N,1))]
    collision, info = env.CheckCollisionRays(rays,target)
    # make sure all normals are the correct sign: pointing outward from the object)
    newinfo = info[collision,:]
    if len(newinfo) > 0:
          newinfo[sum(rays[collision,3:6]*newinfo[:,3:6],1)>0,3:6] *= -1
          approachrays = r_[approachrays,newinfo]
    if normalanglerange > 0:
         theta,pfi = SpaceSamplerExtra().sampleS2(angledelta=directiondelta)
         dirs = c_[cos(theta),sin(theta)*cos(pfi),sin(theta)*sin(pfi)]
         dirs = array([dir for dir in dirs if arccos(dir[2])<=normalanglerange]) # find all dirs within normalanglerange
         if len(dirs) == 0:
                 dirs = array([[0,0,1]])
         newapproachrays = zeros((0,6))
         for approachray in approachrays:
             R = rotationMatrixFromQuat(quatRotateDirection(array((0,0,1)),approachray[3:6]))
             newapproachrays = r_[newapproachrays,c_[tile(approachray[0:3],(len(dirs),1)),dot(dirs,transpose(R))]]
         approachrays = newapproachrays
    return approachrays

def NeighborPoints(point,anglerange,distance):
    return
    
def PointsForCoating(points, distance):
    for point in points:
        point[0:3] = point[0:3]+distance*point[3:6]
    return points

def hat(vec):
    hvec = [[0, -vec[2], vec[1]],[vec[2],0,-vec[0]],[-vec[1], vec[0],0]]
    return hvec

def RotationCalc(a,b): # de a para b
    v = cross(a,b)
    vhat = hat(v)
    cosab = dot(a,b)
    sinab = dot(v,v)
    R = eye(3)+vhat+numpy.matrix_power(vhat,2)*(1-cosab)/(sinab**2)
    R = transpose(R)
    return R

def PointTransform(initialorientation, pointnormal):
    normal = pointnormal[3:6]
    point = pointnormal[0:3]
    R = RotationCalc(initialorientation,normal)
    T = [[R[0][0],R[0][1],R[0][2],point[0]],[R[1][0],R[1][1],R[1][2],point[1]],
         [R[2][0],R[2][1],R[2][2],point[2]],[0,0,0,1]]
    return T

def genVec(normal):
    if dot(normal,[1,0,0])==1:
        v2 = [0,1,0]
    else:
        v2 = [1,0,0]
    v2 = cross(normal,v2)
    v3 = cross(normal,v2)
    return v2,v3

def genTransform(pointnormal):
    normal = pointnormal[3:6]
    point = pointnormal[0:3]
    v1 = normal
    v2,v3 = genVec(normal)
    T = [[v1[0],v2[0],v3[0],point[0]],[v1[1],v2[1],v3[1],point[1]],
         [v1[2],v2[2],v3[2],point[2]],[0,0,0,1]]
    return T
def pose(pointnormal, distance):
    normal = pointnormal[3:6]
    point = pointnormal[0:3]
    pM = point+distance*normal

    v1 = normal
    v1[1]=0
    v1=v1/sqrt(dot(v1,v1))
    v2 = [0,1,0]
    v3 = cross(v1,v2)

    T = [[v1[0],v2[0],v3[0],pM[0]],[v1[1],v2[1],v3[1],pM[1]],
         [v1[2],v2[2],v3[2],pM[2]],[0,0,0,1]]
    
    return T

def WorkspaceOnPose(robotpose, distance, bladepoints):
    Tn = pose(robotpose,distance)
    robot.SetTransform(Tn)
    initial_angles = [ 0,  0,  0,  0,  0,  0]
    robot.SetDOFValues(initial_angles,ikmodel.manip.GetArmIndices())
    reachableRays=zeros((0,6))
    for ray in bladepoints:
        T = genTransform(ray)
        iksol = ikmodel.manip.FindIKSolution(T,IkFilterOptions.CheckEnvCollisions)
        if iksol is not None:
            reachableRays = vstack((reachableRays,ray))
    return reachableRays

def BestBaseDistance(robotpose,initialdistance,bladepoints):
    distance = initialdistance
    reachableRays=zeros((0,6))
    bestDistance = 0
    for i in range(0,181):
        distance = initialdistance + 1.0*i/100
        tempReachableRays = WorkspaceOnPose(robotpose, distance, bladepoints)
        if len(tempReachableRays) > len(reachableRays):
            reachableRays = tempReachableRays
            bestDistance = distance
    Tn = pose(robotpose,bestDistance)
    robot.SetTransform(Tn)        
    return  reachableRays, bestDistance       
    
#MAIN
with env:

    approachrays=PointsToReach(env, target, delta, normalanglerange,directiondelta)
    approachrays=PointsForCoating(approachrays,coatingdistance)
    N = approachrays.shape[0]
    Ttarget = target.GetTransform()

    I=0
    M=N-I     
    #PLOT BLADE POINTS FOR COATING
    gapproachrays = c_[dot(approachrays[I:N,0:3],transpose(Ttarget[0:3,0:3]))+tile(Ttarget[0:3,3],(M,1)),dot(approachrays[I:N,3:6],transpose(Ttarget[0:3,0:3]))]
    approachgraphs = env.plot3(points=gapproachrays[:,0:3],pointsize=5,colors=array((1,0,0)))
    
    #reachableRays = WorkspaceOnPose(pN, robottobladedistance, approachrays)
    #reachableRays, bestDistance = BestBaseDistance(pN, initialdistance, approachrays)
    #print bestDistance
    
    #PLOT REACHABLE POINT
    #grays = c_[dot(reachableRays[:,0:3],transpose(Ttarget[0:3,0:3]))+tile(Ttarget[0:3,3],(len(reachableRays),1)),dot(reachableRays[:,3:6],transpose(Ttarget[0:3,0:3]))]
    #raygraphs = env.plot3(points=reachableRays[:,0:3],pointsize=5,colors=array((0,0,0)))
