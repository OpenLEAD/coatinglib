# Function definition
from numpy import *
from openravepy import *
from openravepy.misc import SpaceSamplerExtra
import time
import math

def CoordString(x): # Chooses the plane of the blade to discretize
    return {
        'b': 0,
        'r': 1,
        't':2,
        'l':3,
        's': 4,
        's2': 5,
    }[x]


def PointsToReach(env, target, delta, normalanglerange,directiondelta,lrtb): # Discretizes the blade in points to coat
    cc = RaveCreateCollisionChecker(env,'ode')
    if cc is not None:
            ccold = env.GetCollisionChecker()
            env.SetCollisionChecker(cc)
            cc = ccold

    target.SetTransform(eye(4))
    ab = target.ComputeAABB()
    p = ab.pos()
    e = ab.extents()#+0.02 # increase since origin of ray should be outside of object
    sides = array(((0,0,e[2],0,0,-1,e[0],0,0,0,e[1],0),
                 (0,0,-e[2],0,0,1,e[0],0,0,0,e[1],0),
                 (0,e[1],0,0,-1,0,e[0],0,0,0,0,e[2]),
                 (0,-e[1],0,0,1,0,e[0],0,0,0,0,e[2]),
                 (e[0],0,0,-1,0,0,0,e[1],0,0,0,e[2]),
                 (-e[0],0,0,1,0,0,0,e[1],0,0,0,e[2])))
    maxlen = 2*sqrt(sum(e**2))+0.03
    approachrays = zeros((0,6))
    x=CoordString(lrtb)
    side=sides[x];
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

    
def PointsForCoating(points, distance): # Points are shifted distance from the blade (to meet requirement of the coating distance)
    for point in points:
        point[0:3] = point[0:3]+distance*point[3:6]
    return points

def hat(vec): # Skew-symmetric matrix
    hvec = array([[0, -vec[2], vec[1]],[vec[2],0,-vec[0]],[-vec[1], vec[0],0]])
    return hvec

def TensorDot(a,b): # a*b.transpose()
     R = array([[a[0]*b[0],a[0]*b[1],a[0]*b[2]],
                      [a[1]*b[0],a[1]*b[1],a[1]*b[2]],
                      [a[2]*b[0],a[2]*b[1],a[2]*b[2]]])
     return R

def RotationCalc(a,b): # Matrix rotation from 'a' to 'b'
    a = a/sqrt(dot(a,a))
    b = b/sqrt(dot(b,b))

    v = cross(a,b)
    sinab = sqrt(dot(v,v))
    vhat = hat(v)
    cosab = dot(a,b)

    if sinab == 0 and cosab == 1:
        R = eye(3)
    elif sinab == 0 and cosab == -1:
        R = -eye(3)
    else:    
        R = eye(3)+vhat+dot(vhat,vhat)*(1-cosab)/(sinab**2)
    return R

def RotationCalc2(a,b): # Matrix rotation on 'a' axis, angle 'b'
     R = eye(3)*cos(b) + hat(a)*sin(b) + (1-cos(b))*TensorDot(a,a)
     return R

def PointTransform(initialorientation, pointnormal):
    normal = pointnormal[3:6]
    point = pointnormal[0:3]
    R = RotationCalc(initialorientation,normal)
    T = [[R[0][0],R[0][1],R[0][2],point[0]],[R[1][0],R[1][1],R[1][2],point[1]],
         [R[2][0],R[2][1],R[2][2],point[2]],[0,0,0,1]]
    return T

def genVec(normal): # Given a direction (normal), generates two arbitrary vectors to compose a space
    if dot(normal,[1,0,0])==1:
        v1 = [0,1,0]
    else:
        v1 = [1,0,0]
    v1 = cross(normal,v1)
    v1 = v1/sqrt(dot(v1,v1))
    v3 = cross(normal,v1)
    v3 = v3/sqrt(dot(v3,v3))
    return v1,v3

def genTransform(pointnormal, facevector): # Generates the matrix T, given point, normal vector and vector to face normal
    normal = pointnormal[3:6]
    point = pointnormal[0:3]

    R = RotationCalc(facevector,normal)

    T = [[R[0][0],R[0][1],R[0][2],point[0]],[R[1][0],R[1][1],R[1][2],point[1]],
         [R[2][0],R[2][1],R[2][2],point[2]],[0,0,0,1]]
    return T

def normalmean(bladepoints): # Computes the average of all normal vectors of the blade
    normal = [0,0,0]
    for ray in bladepoints:
        normal = [normal[0]+ray[3],normal[1]+ray[4],
                  normal[2]+ray[5]]
    normal = array(normal)
    divider = sqrt(dot(normal,normal))
    mean = normal/divider
    return mean

def pose(pointnormal, distance, bladepoints): # Computes robot position in respect of average of all normal vectors of the blade
    #normal = pointnormal[3:6]
    normal = normalmean(bladepoints)
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

def poseDummy(pointnormal, distance): # Computes robot position in respect of given normal vector, point and distance
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

def Rz(T,theta): # Rotate robot in Z axis, given matrix T of the robot and angle theta to rotate
    R = [[T[0][0],T[0][1],T[0][2]],[T[1][0],T[1][1],T[1][2]],
         [T[2][0],T[2][1],T[2][2]]]
    v1 = array([cos(theta), -sin(theta),0])
    v1R = dot(v1,R)
    v2 = array([sin(theta),cos(theta),0])
    v2R = dot(v2,R)
    v3 = array([0,0,1])
    v3R = dot(v3,R)

    T = [[v1R[0],v1R[1],v1R[2],T[0][3]],[v2R[0],v2R[1],v2R[2],T[1][3]],
         [v3R[0],v3R[1],v3R[2],T[2][3]],[0,0,0,1]]
    return T

def RzM(T,theta): # Rotate robot in Z axis (world frame), given matrix T of the robot and angle theta to rotate
    R = [[T[0][0],T[0][1],T[0][2]],[T[1][0],T[1][1],T[1][2]],
         [T[2][0],T[2][1],T[2][2]]]
    R = array(R)
    v1 = array([cos(theta), sin(theta),0])
    v1R = dot(v1,R.T)
    v2 = array([-sin(theta),cos(theta),0])
    v2R = dot(v2,R.T)
    v3 = array([0,0,1])
    v3R = dot(v3,R.T)

    T = [[v1R[0],v2R[0],v3R[0],T[0][3]],[v1R[1],v2R[1],v3R[1],T[1][3]],
         [v1R[2],v2R[2],v3R[2],T[2][3]],[0,0,0,1]]
    return T

def RxM(T,theta):  # Rotate robot in X axis (world frame), given matrix T of the robot and angle theta to rotate
    R = [[T[0][0],T[0][1],T[0][2]],[T[1][0],T[1][1],T[1][2]],
         [T[2][0],T[2][1],T[2][2]]]
    R = array(R)
    v1 = array([1, 0, 0])
    v1R = dot(v1,R.T)
    v2 = array([0, cos(theta), sin(theta)])
    v2R = dot(v2,R.T)
    v3 = array([0, -sin(theta), cos(theta)])
    v3R = dot(v3,R.T)

    T = [[v1R[0],v2R[0],v3R[0],T[0][3]],[v1R[1],v2R[1],v3R[1],T[1][3]],
         [v1R[2],v2R[2],v3R[2],T[2][3]],[0,0,0,1]]
    return T


def defaultSolve(ikmodel, T):
    iksol = ikmodel.manip.FindIKSolutions(T,IkFilterOptions.CheckEnvCollisions)
    if len(iksol)>0:
        #iksolpositive = []
        #for sol in iksol:
           #for i in range(0,len(sol)):
                #if sol[i]<0:sol[i]+=2*math.pi
           #iksolpositive.append(sol)     
        return iksol, True
    else:
        return iksol, False
    
def KinematicSolve(bladepoint,ikmodel,facevector,coatingdistancetolerance=0): # Solve inverse kinematics for specific point, normal vector, robot, and vector to face normal vector of the blade
    T = genTransform(bladepoint,facevector)
    iksol, answer = defaultSolve(ikmodel, T)
    if coatingdistancetolerance!=0:
        bladepoint[0:3] = bladepoint[0:3]+coatingdistancetolerance*bladepoint[3:6]
        T = genTransform(bladepoint,facevector)
        iksol2, answer2 = defaultSolve(ikmodel, T)
    iksol = concatenate((iksol,iksol2))
    answer = answer or answer2
    return iksol, answer        

def AllKinematicSolve(bladepoints,ikmodel,facevector,coatingdistancetolerance=0): # Solve inverse kinematics for all specific points, normal vectors, robot, and vector to face normal vector of the blade
     reachableRays=zeros((0,6))
     iksolList = []
     indexlist = zeros((len(bladepoints),1),dtype=bool)
     i=0
     for ray in bladepoints:
        iksol, index = KinematicSolve(ray,ikmodel,facevector,coatingdistancetolerance)
        if len(iksol)>0:
            reachableRays = vstack((reachableRays,ray))
            for ik in iksol:
                 iksolList.append(ik)
        indexlist[i] = index    
        i+=1    
     return reachableRays, iksolList, indexlist

def WorkspaceOnPose(robotpose, distance, bladepoints,robot,ikmodel,facevector,theta,coatingdistancetolerance=0): # Pose the robot and solve inverse kinematics for all specific points, normal vectors, robot, and vector to face normal vector of the blade.
    thetaX = theta[0]
    thetaY = theta[1]
    thetaZ = theta[2]
    #Tn = pose(robotpose,distance,bladepoints)
    Tn = poseDummy(robotpose,distance)
    Tn = RzM(Tn,thetaZ)
    Tn = RxM(Tn,thetaX)
    robot.SetTransform(Tn)
    reachableRays=zeros((0,6))
    iksolList = []
    indexlist = zeros((len(bladepoints),1),dtype=bool)
    i=0
    for bladepoint in bladepoints:
          ik, index = KinematicSolve(bladepoint,ikmodel,facevector,coatingdistancetolerance)
          if index:
               iksolList.append(ik)
               reachableRays = vstack((reachableRays,bladepoint))
          indexlist[i]=index
          i+=1
    return reachableRays, iksolList, indexlist

def BestBaseDistance(robotpose,initialdistance,bladepoints,robot,ikmodel,facevector,theta,coatingdistancetolerance=0): # Computes best robot position to coat blade points
    distance = initialdistance
    reachableRays=zeros((0,6))
    bestDistance = 0
    for i in range(0,181):
        distance = initialdistance + 1.0*i/100
        tempReachableRays, iksolList = WorkspaceOnPose(robotpose, distance, bladepoints,robot,ikmodel,facevector,theta,coatingdistancetolerance)
        if len(tempReachableRays) > len(reachableRays):
            reachableRays = tempReachableRays
            bestDistance = distance
        print 'iterator:', i, ', number of points:', len(tempReachableRays), ', distance:', distance, ', bestDistance:', bestDistance
    Tn = pose(robotpose,bestDistance,bladepoints)
    robot.SetTransform(Tn)        
    return  reachableRays, bestDistance

def robotPath(iksolList, timesleep,robot,ikmodel): # Show robot motion
    for sol in iksolList:
        robot.SetDOFValues(sol,ikmodel.manip.GetArmIndices())
        time.sleep(timesleep)
    return 

def NotCoatedRays(approachrays,reachableRays): # Computes points that were not coated
     s1 = set(map(tuple,approachrays))
     s2 = set(map(tuple,reachableRays))
     s3 = s1-s2
     notReachableRays = array(map(tuple,s3))
     return notReachableRays

def genRays(numberofangles,normal,p,newpoint,distance): # Generate points that are equivalent to reference point due to tolerance
     newbladepoints=zeros((0,6))
     k = 1.0*2*pi/numberofangles
     for i in range(0,numberofangles):
          alfa = k*i
          Rv1alfa = RotationCalc2(normal,alfa)
          tempP = dot(p,transpose(Rv1alfa))

          point = newpoint+distance*tempP
          newnormal = tempP
          ray = concatenate((point,newnormal))
          
          newbladepoints = vstack((newbladepoints,ray))
     return newbladepoints

def TryCoatnonNormalCoatPoint(pointnormal, distance, numberofangles, tolerance, ikmodel, facevector,coatingdistancetolerance=0.01): # try to coat point that were not coated with tolerance
     normal = pointnormal[3:6]
     point = pointnormal[0:3]
     newpoint = point-distance*normal
     tolerance = 1.0*pi*tolerance/180

     v1 = normal 
     v2, v3 = genVec(v1)

     Rv3tol = RotationCalc2(v3,tolerance)
     p = dot(v1,transpose(Rv3tol))

     newbladepoints=genRays(numberofangles,v1,p,newpoint,distance)
     
     reachableRays, iksolList, _ = AllKinematicSolve(newbladepoints,ikmodel,facevector,coatingdistancetolerance)

     if len(reachableRays)>0:
          return True, iksolList
     else:
          return False, iksolList

def AllExtraCoating2(approachrays,indexreachableRays,distance, numberofangles, tolerance, ikmodel, facevector,coatingdistancetolerance=0.01): # try to coat all points with tolerance
     AllreachableRays=zeros((0,6))
     AlliksolList = []
     indexlist = zeros((len(approachrays),1),dtype=bool)
     i=0
     for index in indexreachableRays:
          if not index:
               answer, iksol = TryCoatnonNormalCoatPoint(approachrays[i], distance, numberofangles, tolerance, ikmodel, facevector,coatingdistancetolerance)
               if answer:
                    AllreachableRays = vstack((AllreachableRays,approachrays[i]))
                    AlliksolList.append(iksol)
                    indexlist[i]=1
               else:
                    indexlist[i]=0
          else:
               indexlist[i]=0
          i+=1     
     return AllreachableRays, AlliksolList, indexlist      

def IndexToPoints(bladepoints,indexlist): # converts boolean array of [coated,non coated] blade points to [x,y,z,normal]
     reachableRays=zeros((0,6))
     i=0
     for index in indexlist:
          if index:
               reachableRays = vstack((reachableRays,bladepoints[i]))
          i+=1
     return reachableRays     

def clearPathIntersection(fullindexlist): # clears intersections of coating process
     for i in range(0,len(fullindexlist)):
          for j in range(i+1,len(fullindexlist)):
               fullindexlist[i]=fullindexlist[i]&~fullindexlist[j]
     return fullindexlist           

def iksolsSort(indexlist1,indexlist2,iksollist1,iksollist2,ikmodel): # concatenates coating solutions without and with tolerance
     iksollist = []
     index1=0
     index2=0
     for i in range(0,len(indexlist1)):
          if indexlist1[i]:
               iksollist.append(iksollist1[index1])
               index1+=1
          elif indexlist2[i]:
               iksollist.append(iksollist2[index2])
               index2+=1
     return iksollist          

def nearPointsByDistance(arraylist): # calculates nearest points by distance
     n = len(arraylist)
     nearlist = []
     for i in range(0,n):
          nearpoints = []
          nearpoints.append(arraylist[i])
          for array in arraylist:
               arrayp = array([array[0],array[1],array[2]])
               arraylistp = array([arraylist[i][0],arraylist[i][1],arraylist[i][2]])
               dist = linalg.norm(arrayp-arraylistp)
               if dist<0.1 and dist!=0:
                    nearpoints.append(array)
          nearlist.append(nearpoints)
          #print str(i)+'/'+str(n)
     return nearlist

def nearPointsByNumberOfPoints(arraylist): # find the 8 nearest points 
     NumberOfPoints = 8
     nearlist = []
     n=len(arraylist)
     k=0
     for arrays in arraylist:
          nearpoints = []
          nearpoints.append(arrays)
          dist = []
          for i in range(0,len(arraylist)):
               arrayp = array([arrays[0],arrays[1],arrays[2]])
               arraylistp = array([arraylist[i][0],arraylist[i][1],arraylist[i][2]])
               dist.append(linalg.norm(arrayp-arraylistp))
          thelist = [x for (y,x) in sorted(zip(dist,arraylist), key=lambda pair: pair[0])]     
          for j in range(1,NumberOfPoints+1):
               nearpoints.append(thelist[j])
          nearlist.append(nearpoints)
          k+=1
          print str(k)+'/'+str(n)
     return nearlist
     
def projectPointInPlane(pointtoproject,pointinplane,normalofplane): # projects points in given plane
     v = pointtoproject-pointinplane
     d = dot(v,normalofplane)
     pointprojection = pointtoproject - d*normalofplane
     return pointprojection

def angularDistance(referencepoint,pointprojections,normalplane): # calculates angular difference between first vector and others
     angles=[]
     v = pointprojections[0]-referencepoint
     vnorm = sqrt(dot(v,v))
     for pointprojection in pointprojections:
          v2 = pointprojection-referencepoint
          thecross = cross(v,v2)
          v2norm = sqrt(dot(v2,v2))
          thecos = 1.0*dot(v,v2)/(vnorm*v2norm)
          thesin = 1.0*dot(thecross/(vnorm*v2norm),normalplane)
          angles.append(math.atan2(thesin,thecos))
     return angles

def computeAngularDistances(points): 
     reference = points[0]
     pointprojections = []
     referencepoint = array([reference[0],reference[1],reference[2]])
     normalplane = array([reference[3],reference[4],reference[5]])
     n = len(points)
     for i in range(1,n):
          pointnormal = points[i]
          point = array([pointnormal[0],pointnormal[1],pointnormal[2]])
          pointprojections.append(array(projectPointInPlane(point,referencepoint,normalplane)))
     angles = angularDistance(referencepoint,pointprojections,normalplane)
     return angles#, pointprojections

def computeAllAngularDistances(allpoints):
     allangles = []
     n=len(allpoints)
     i=0
     for points in allpoints:
          allangles.append(computeAngularDistances(points))
          i+=1
          #print str(i)+'/'+str(n)
     return allangles     

def clusteringNearestAngles(angles,nearlist):
     reference = nearlist[0]
     nearlist = delete(nearlist,0,0)
     sortednearlist = [x for (y,x) in sorted(zip(angles,nearlist), key=lambda pair: pair[0])]
     angles = sorted(angles)
     clusterlistangles = []
     clusterlistpoints = []
     initialangle = angles[0]
     nearangles = []
     nearpoints = []
     for i in range(0,len(angles)):
          delta = abs(initialangle-angles[i])
          if delta<math.pi/4 or abs(delta-2*math.pi)<math.pi/4:
               nearangles.append(angles[i])
               nearpoints.append(sortednearlist[i])
          else:
               nearpoints.append(reference)
               clusterlistangles.append(nearangles)
               clusterlistpoints.append(nearpoints)
               nearangles = []
               nearpoints = []
               nearangles.append(angles[i])
               nearpoints.append(sortednearlist[i])
          initialangle = angles[i]     
     return clusterlistangles,clusterlistpoints
               
def clusteringAllNearestAngles(allangles,allnearlist):
     allclusterlistangles=[]
     allclusterlistpoints=[]
     n = len(allnearlist)
     for i in range(0,len(allnearlist)):
          clusterlistangles,clusterlistpoints = clusteringNearestAngles(allangles[i],allnearlist[i])
          allclusterlistangles.append(clusterlistangles)
          allclusterlistpoints.append(clusterlistpoints)
          #print str(i)+'/'+str(n)
     return allclusterlistangles,allclusterlistpoints     
          
def pairingAngles(angles,points):
     couplesangles=[]
     triopoints=[]
     reference = points[0]
     points = delete(points,0,0)
     while(len(angles)):
          oneangle = angles[0]
          onepoint = points[0]
          angles = delete(angles,0,0)
          points = delete(points,0,0)
          for j in range(0,len(angles)):
               if abs(oneangle-angles[j])-math.pi<0.05:
                    couplesangles.append(array([oneangle,angles[j]]))
                    triopoints.append(array([reference,onepoint,points[j]]))
                    angles = delete(angles,j,0)
                    points = delete(points,j,0)
                    break
          if j==len(angles)-1: return False, couplesangles, triopoints
     return True, couplesangles, triopoints

def pairingAllAngles(allangles,allpoints):
     allcouplesangles=[]
     alltriopoints=[]
     for i in range(0,len(allangles)):
          test, couplesangles, triopoints = pairingAngles(allangles[i],allpoints[i])
          if test:
               allcouplesangles.append(couplesangles)
               alltriopoints.append(triopoints)
     return allcouplesangles, alltriopoints
               
def deltaTCalc(alltriopoints):
     deltasT = []
     v = 1.0*40/60 #m/s
     for triopoints in alltriopoints:
          deltaT = []
          for triopoint in triopoints:
               r = array([triopoint[0][0],triopoint[0][1],triopoint[0][2]])
               p1 = array([triopoint[1][0],triopoint[1][1],triopoint[1][2]])
               p2 = array([triopoint[2][0],triopoint[2][1],triopoint[2][2]])
               d1=linalg.norm(r-p1)
               d2=linalg.norm(r-p2)
               deltaT.append(array([d1/v,d2/v]))
          deltasT.append(deltaT)
     return deltasT

def minimize2AngularVelocity(iksols1,iksols2):
     minis = []
     minimalw = 700
     for iksol1 in iksols1:
          for iksol2 in iksols2:
               tempmin = min(abs(iksol1-iksol2))
               if tempmin<minimalw:
                    minimalw = tempmin
                    minis = [iksol1,iksol2]
     return minis

def minimize1AngularVelocity(iksol1,iksols2):
     minis = []
     minimalw = 700
     for iksol2 in iksols2:
          tempmin = min(abs(iksol1-iksol2))
          if tempmin<minimalw:
               minimalw = tempmin
               minis = [iksol1,iksol2]
     return minis

def alphaCalc(triopoint,pN, robottobladedistance,robot,ikmodel,facevector,theta,numberofangles,tolerance,coatingdistance,deltaT):
     _, thetas1, indexlist = WorkspaceOnPose(pN, robottobladedistance, triopoint,robot,ikmodel,facevector,theta)
     _, thetas2, indexlist2  = AllExtraCoating2(triopoint,indexlist,coatingdistance,numberofangles,tolerance,ikmodel,facevector)
     thetas=iksolsSort(indexlist,indexlist2,thetas1,thetas2,ikmodel)
     minis=minimize2AngularVelocity(thetas[0],thetas[1])
     minis2=minimize1AngularVelocity(minis[0],thetas[2])
     Thetas = array([minis[1],minis[0],minis2[1]])
     dtheta1 = minis[1]-minis[0]
     dtheta2 = minis2[0]-minis2[1]
     omega1 = dtheta1/deltaT[0]
     omega2 = dtheta2/deltaT[1]
     alpha = (omega1-omega2)/((deltaT[0]+deltaT[1])/2)
     return alpha,omega1,omega2,Thetas

def AllalphaCalc(alltriopoints,pN, robottobladedistance,robot,ikmodel,facevector,theta,numberofangles,tolerance,coatingdistance,deltasT):               
     alphas = []
     omegas = []
     Thetas = []
     n = len(alltriopoints)
     i=0
     for triopoints in alltriopoints:
         j=0
         alpha = []
         omega = []
         Theta = []
         for triopoint in triopoints:
             alpha0,omega1,omega2,thetas = alphaCalc(triopoint,pN, robottobladedistance,robot,ikmodel,facevector,theta,numberofangles,tolerance,coatingdistance,deltasT[i][j])
             alpha.append(alpha0)
             omega.append(array([omega1,omega2]))
             Theta.append(thetas)
             j+=1
         Thetas.append(Theta)    
         alphas.append(alpha)
         omegas.append(omega)
         i+=1
         print str(i)+'/'+str(n)
     return alphas, omegas,Thetas
 
def calculateOmegasbyJacobian(robot,ikmodel,manip,thetas,alltriopoints,deltasT):
     NewOmegas=[]
     for i in range(0,len(thetas)):
         NewOmega = []
         for j in range(0,len(thetas[i])):
             Omega = []
             robot.SetDOFValues(thetas[i][j][1],ikmodel.manip.GetArmIndices())
             Jpos = manip.CalculateJacobian()
             Jori = manip.CalculateAngularVelocityJacobian()
             J = concatenate((Jpos,Jori))
             invJ = linalg.pinv(J)
             r = array([alltriopoints[i][j][0][0],alltriopoints[i][j][0][1],alltriopoints[i][j][0][2]])
             p1 = array([alltriopoints[i][j][1][0],alltriopoints[i][j][1][1],alltriopoints[i][j][1][2]])
             p2 = array([alltriopoints[i][j][2][0],alltriopoints[i][j][2][1],alltriopoints[i][j][2][2]])
             v1 = (p1-r)/deltasT[i][j][0]
             v2 = (r-p2)/deltasT[i][j][1]
             v1=concatenate((v1,array([0,0,0])))
             v2=concatenate((v2,array([0,0,0])))
             Omega.append(dot(invJ,v1))
             Omega.append(dot(invJ,v2))
             NewOmega.append(Omega)
         NewOmegas.append(NewOmega)     
     return NewOmegas
