from numpy import array, dot, cross, outer, eye, sum, sqrt, exp
from numpy import random, transpose, zeros, linalg, multiply, inf
from numpy import ndindex, linspace, power,  ceil, floor, einsum
from numpy import argsort, argmin, argmax, linspace, power, arange
from numpy import ones, maximum, minimum, round, sign, vander, unravel_index
from numpy import cos as npcos
from numpy import sin as npsin
from numpy import tan as nptan
from math import cos, sin, pi, isnan, acos, atan2
from abc import ABCMeta, abstractmethod
from copy import copy, deepcopy
from collections import deque
from scipy.spatial import KDTree
from scipy.ndimage.interpolation import shift
from itertools import cycle, islice
from scipy.linalg import logm, expm
from openravepy import quatFromRotationMatrix, matrixFromQuat

_base_module_precision = 0.2
_p_step = 0.1

class KinBodyError(Exception):    
    def __init__(self):
        Exception.__init__(self, "object is not a KinBody.")

def csv_compute_hull(csv_file, bases):
    import csv
    from rail_place import RailPlace
    from scipy.spatial import ConvexHull
    positions = []
    angles = []
    
    with open(csv_file, 'rb') as csvfile:
	spamreader = csv.reader(csvfile)
	titles = spamreader.next()
	rotor_col = titles.index('Rotor Angle')
	base_col = titles.index('Base Position')
	for row in spamreader:
            positions += [int(row[base_col])]
            angles += [float(row[rotor_col])]
    
    regions = {}
    for position,angle in zip(positions,angles):
	regions[angle] = regions.get(angle,[]) + [position]
                       
    onipositions = reduce(lambda x,y: set(x) & set(y), regions.values())
    onixy = [ RailPlace(bases[psa]).getXYZ()[0:2] for psa in onipositions ]
    hull2D = ConvexHull(onixy)
    points2D = array(onixy)

    return points2D, hull2D

def base_points_by_angle( points2D, angle, turb ):
    from rail_place import xya2ps, RailPlace
    from itertools import repeat
    step = _base_module_precision * cos(angle)
    y = round(points2D[:,1]/step) * step
    p,s = xya2ps(points2D[:,0], y, angle)

    viable = []
    for psa in zip(p,s,repeat(angle)):
        rp = RailPlace(psa)
        turb.place_rail(rp)
        turb.place_robot(rp)
        viable += [not( turb.check_rail_collision() &
                   turb.check_robotbase_collision())]


    return viable
    
def plot_hull(points2D, hull2D, plt = None):
    if plt == None:
        import matplotlib.pyplot as plt
        
    rndhull = list(hull2D.vertices)+ [hull2D.vertices[0]]

    plt.scatter(array(points2D)[:,0],array(points2D)[:,1])
    plt.plot(array(points2D)[rndhull,0], array(points2D)[rndhull,1], 'r--', lw=2)
    plt.plot(array(points2D)[rndhull[:-1],0], array(points2D)[rndhull[:-1],1], 'ro')
    return plt

def secondary_positions_by_angle(points2D, angle):
    step = _base_module_precision * cos(angle)
    points = array(points2D)
    near = min(points[:,1])
    far = max(points[:,1])

    p0 = ceil(near/step)
    pf = floor(far/step)
    return max(int(pf-p0),0)
    
    
def base_region_by_angle( polygon_vertex2D, angle, turb = None ):
    
    pv = array(polygon_vertex2D)
    
    if all(pv[-1] == pv[0]):
        pv = pv[:-1]
        
    step = _base_module_precision * cos(angle)
    near = min(pv[:,1])
    far = max(pv[:,1])

    p0 = ceil(near/step)
    pf = floor(far/step)

    H = array( [ array([[1,0],pv[n]-pv[n-1]]).T for n in range(len(pv)) ] )

    border = []
    
    for sy in arange(p0*step, (pf+1)*step, step):
        O = array([0,sy])
        dt = linalg.solve(H, pv - O)
        hits = dt[(dt[:,1]>=0) & (dt[:,1]<=1)]
        hits = sorted(hits[:,0])

        for segment in zip(hits[0::2],hits[1::2]):
            segborder = []
            first, last = True, False
            for point in range(segment[0],segment[1]+_p_step,_p_step):
                p,s = xya2ps(point, sy, angle)
                rp = RailPlace((p,s,angle))
                turb.place_rail(rp)
                turb.place_robot(rp)
                viable = not( turb.check_rail_collision() &
                           turb.check_robotbase_collision())
                if viable and first:
                    segborder += [(point,sy)]
                    first = False
                #if not#CONTINUE
            
        
def direction_in_halfplane(rays,direction):
    """
    Filter points on the half plane defined by
    normal plane to direction ( dot(ray[3:6],direction) > 0 )
    """
    return array(rays)[ (dot(array(rays)[:,3:6],direction)) > 0 ]
    

def annulus_distribution(N,r,R, origin = None, dim = 2):
    """
    return N samples from a uniform distribution over a
    annulus with inner radius r and outer radius R.

    """
    if origin is None:
        origin = zeros((1,dim))
    
    direction = random.randn(N,dim)
        
    rho = random.rand(N)
    rho = power((R**dim - r**dim)*rho + r**dim,1./dim)
    
    samples = einsum('ij,i->ij',direction,rho/linalg.norm(direction, axis=1))

    return samples + origin


def fast_poisson_disk(r, limits, k = 30, points = None):
    n = len(limits)
    cellfreq = sqrt(n)/r
    delta = dot(limits,[-1, 1])
    gridsize = ceil(delta*cellfreq)
    grid = - ones(gridsize.astype('int')).astype('int')
    if points is None:
        x0 = random.rand(n)
        grid[tuple(floor(x0*gridsize).astype('int'))] = 0
        x = [limits[:,0]+delta*x0]
        activelist = [0]
    else:
        x = []
        activelist = []
        for point in points:
            grid_place = floor((point-limits[:,0])*gridsize/delta).astype('int')
            if grid[tuple(grid_place)] == -1:
                grid[tuple(grid_place)] = len(x)
                activelist += [len(x)]
                x += [point]

    while len(activelist)>0:
        i = activelist[random.randint(len(activelist))]
        samples = annulus_distribution(k,r,2*r,x[i],n)
        placed = False
        for sample in samples:
            grid_place = floor((sample-limits[:,0])*gridsize/delta).astype('int')
            if any(grid_place < 0) or any(grid_place >= grid.shape) or (grid[tuple(grid_place)] != -1):
                continue

            for idx in ndindex(tuple([3]*n)):
                try:
                    point_index = grid[tuple(array(idx)-1+grid_place)]
                except IndexError:
                    continue
                if point_index == -1:
                    continue
                if linalg.norm(sample - x[point_index]) < r:
                    break
            else:
                grid[tuple(grid_place)] = len(x)
                activelist += [len(x)]
                x += [sample]
                placed = True
        if not placed:
            del activelist[activelist.index(i)]
    return x

def update_rail_region(cfg, db_bases_to_num, db_visited_bases, distance = None, density = None, extra_points = None):
    """
    cfg = TurbineConfig with new intended limits (should be a higher limit)
    db_bases_to_num, db_visited_bases = the data bases to be updated (std dictionary format)
    
    distance, density, extra_points = Information about sampled data points in order of priority
                                    (if one info is provided the others are ignored)

    distance = expected distance between base points
    density  = expected density of base point (points per area squared)
    extra_points = number of desired extra points
    If None is given, it tries to estimate the actual distance between points and use it.
    
    """
    from time import time
    import cPickle
    from scipy.spatial.distance import pdist
    from rail_place import RailPlace, _rand_angle

    points = [RailPlace(base).getXYZ()[0:2] for base in db_bases_to_num.keys()]
    limits = array([[cfg.environment.x_min, cfg.environment.x_max],
                   [cfg.environment.y_min, cfg.environment.y_max]])
    
    if distance is None:
        if density is None:
            if extra_points is None:
                dists = pdist(points)
                def sqr2cond(i, j, n): # i != j - Square to Condensated Matrix
                    if i < j:
                        i, j = j, i
                    return n*j - j*(j+1)/2 + i - 1 - j

                dist_cum = 0
                N = len(points)
                
                for n,item in enumerate(extra_points):
                    inx = map(lambda x: sqr2cond(x,n,N), range(n)+range(n+1,N))
                    dist_cum += min(dists[inx])

                density = (N*1./dist_cum)**2/sqrt(2)

            else:
                delta_x, delta_y = dot(limits, [-1,1])
                density = (len(points)+extra_points)*1./(delta_x*delta_y)
                
        distance = sqrt(1./(density*sqrt(2)))

    
    full_points = fast_poisson_disk(distance,limits,30,points)

    alpha_min = cfg.environment.rail_angle_mean - cfg.environment.rail_angle_limit
    alpha_max = cfg.environment.rail_angle_mean + cfg.environment.rail_angle_limit


    x,y = transpose(full_points)
    y_max = (cfg.environment.x_max - x)/nptan(alpha_min)
    y_min = (cfg.environment.x_min - x)/nptan(alpha_min)
    y_max = minimum(y_max,cfg.environment.y_max)
    y_min = maximum(y_min,cfg.environment.y_min)
    x = x[(y<y_max) & (y>y_min)]
    y = y[(y<y_max) & (y>y_min)]


    spoints = set([tuple(points) for points in points])
    sfullpoints = set([tuple(xpoints) for xpoints in zip(x,y)])

    diffpoints = sfullpoints - spoints
    if len(diffpoints)==0:
        raise ValueError('config file generate no new points.')
    rdiffpoints = array([d for d in diffpoints])

    xdif = transpose(rdiffpoints)[0]
    ydif = transpose(rdiffpoints)[1]

    random.seed(int(time()+int((time()-int(time()))*10000)))
    alpha = random.rand(len(xdif))

    alpha = _rand_angle(cfg, xdif, ydif, alpha)

    S = ydif/npcos(alpha)
    P = xdif + S*npsin(alpha)

    raildif = [RailPlace((p,s,a)) for p,s,a in zip(P,S,alpha)]

    for rail in raildif:
	db_bases_to_num[tuple(rail.getPSAlpha())] = len(db_bases_to_num)

    for rail in raildif:
	db_visited_bases[db_bases_to_num[tuple(rail.getPSAlpha())]] = False

    return db_bases_to_num, db_visited_bases

def central_difference(turbine, joints_trajectory, trajectory_index):

    # 'j' for joints, so it doesnt clumpsy the equations
    j = deque(joints_trajectory)
    j.rotate(-trajectory_index)

    h = turbine.config.model.trajectory_step
    
    # Joints velocity - Central Difference (h**6 order error)
    w = ( (j[3]-j[-3]) + 9*(j[-2]-j[2]) + 45*(j[1]-j[-1]) )/(60.0*h)

    # Joints acceleration - Central Difference (h**6 order error)
    alpha = ( 2*(j[-3]+j[3]) - 27*(j[-2]+j[2]) + 270*(j[-1]+j[1]) - 490*j[0] )/(180.0*h**2)

    return w, alpha

def curvepoint(s1, s2, p0, tol=1e-4):
    """
    Find a point on s1 and s2 (intersection between surfaces) near p0.
    The algorithm is described in the paper:
    'Tracing surface intersections - C.L. Bajaj, C.M. Hoffmann, R.E. Lynch,
    and J.E.H. Hopcroft.'
    Available at:
    https://www.cs.purdue.edu/homes/cmh/distribution/papers/Geometry/geo6.pdf

    Keyword arguments:
    s1 -- surface 1.
    s2 -- surface 2.
    p0 -- initial point (3D - without normal vector).
    tol -- tolerance, stop criteria.
    """
    
    tol2 = tol**2
    while True:
        df1 = s1.df(p0); f1 = s1.f(p0)
        df2 = s2.df(p0); f2 = s2.f(p0)

        crossdf = cross(df1,df2)
        if isnan(crossdf[0]) or isnan(crossdf[1]) or isnan(crossdf[2]):
            raise ValueError('Outside model')
        if (abs(crossdf - array([0,0,0]))<=1e-5).all():
            raise ValueError('Gradients are colinear')
        
        df1df2 = dot(df1,df2); df1df1 = dot(df1,df1); df2df2 = dot(df2,df2)
        beta = (-f1*df1df2+f2*df1df1)/(df1df2**2-df1df1*df2df2)
        alpha = (-f1*df2df2+f2*df1df2)/(df1df1*df2df2-df1df2*df1df2)
        dk = alpha*df1+beta*df2
        p1 = p0+dk
        if dot(dk,dk)<tol2:
            grad = s1.df(p1)
            grad = grad/sqrt(dot(grad,grad))
            p1 = array([p1[0],p1[1],p1[2],
                        grad[0],grad[1],grad[2]])
            return p1
        else:
            p0=p1

def surfaces_tangent(ray, s2):
    """ Find the tangent of two surfaces in point/normal ray 

    Keyword arguments:
    ray -- 6x1 vector, point and normal in surface 1.
    s2 -- surface 2.
    """
    
    tan = cross(ray[3:6],s2.df(ray))
    tan = tan/linalg.norm(tan)
    return tan

def hat(vec):
    """ Skew-symmetric matrix """
    
    hvec = array([[0, -vec[2], vec[1]],
                  [vec[2],0,-vec[0]],
                  [-vec[1], vec[0],0]])
    return hvec

def Rab(a,b):
    """ Matrix rotation from 'a' to 'b' """
    a = a/linalg.norm(a)
    b = b/linalg.norm(b)

    v = cross(a,b)
    sinab = sqrt(dot(v,v))
    vhat = hat(v)
    cosab = dot(a,b)

    if cosab == -1:
        return Raxis(compute_perpendicular_vector(a), pi)
    else:    
        R = eye(3)+vhat+dot(vhat,vhat)*1.0/(1+cosab)
    return R

def Raxis(a,theta):
    """ Matrix rotation on 'a' axis, angle 'theta' """
    
    a = a/linalg.norm(a)    
    R = eye(3)*cos(theta) + hat(a)*sin(theta) + (1-cos(theta))*outer(a, a)
    return R

def rotate_trajectories(trajectories, T):
    """
    Rotate all points in trajectories and return the trajectories.

    Keyword arguments:
    T -- the homogeneous transform matrix 4x4.
    """

    R = T[0:3,0:3]
    for i in range(0,len(trajectories)):
        if len(trajectories[i])>0:
            traj = array(trajectories[i])
            Ra = zeros(traj.shape)
            
            Ra[:,0:3] = dot(traj[:,0:3], transpose(R))
            Ra[:,3:6] = dot(traj[:,3:6], transpose(R))
            
            trajectories[i] = Ra.tolist()
    return trajectories

def rotate_points(trajectories, T=[]):
    """
    Rotate all points in trajectories and return the trajectories.

    Keyword arguments:
    T -- the homogeneous transform matrix 4x4.
    """

    rtrajectories = deepcopy(trajectories)
    if len(T)==0:
        return rtrajectories
    
    R = T[0:3,0:3]
    for i in range(0,len(rtrajectories)):
        traj = [x for x in rtrajectories[i] if len(x)>0]
        if len(traj)==0:
            continue
        traj = array(traj)
        Ra = zeros(traj.shape)
        
        Ra[:,0:3] = dot(traj[:,0:3], transpose(R))

        counter = 0
        for j in range(len(rtrajectories[i])):
            if len(rtrajectories[i][j])>0:
                rtrajectories[i][j] = Ra[counter]
                counter+=1
            else:
                continue
    return rtrajectories
        
def trajectory_verticalization(rays, step = 1):
    counter = 0
    new_rays = []
    while True:
        new_ray = [] 
        for i in range(0,len(rays)):
            if len(rays[i])>counter:
                new_ray.append(rays[i][counter])
        if len(new_ray)>0:
            new_rays.append(new_ray)
            counter+=step
        else: break
    return new_rays

def uneven_trajectory_verticalization(parallels, distance=3e-3):
##    for i in range(len(parallels)):
##        parallels[i] = filter_trajectory(parallels[i], distance)
    mid_parallel = filter_trajectory(parallels[len(parallels)/2], distance)
    verticalized_parallels=zeros((len(mid_parallel),len(parallels),len(mid_parallel[0])))
    verticalized_parallels[:,len(parallels)/2] = mid_parallel
    
    for i in arange(verticalized_parallels.shape[1]/2,0,-1):
        dists = zeros((verticalized_parallels.shape[0],len(parallels[i-1])))
        for j in range(verticalized_parallels.shape[0]):
            dists[j] = linalg.norm(array(parallels[i-1])[:,0:3]-verticalized_parallels[j][i][0:3],axis=1)
        for _ in range(len(dists)):
            j,k = unravel_index(argmin(dists),dists.shape)
            verticalized_parallels[j][i-1] = parallels[i-1][k]
            dists[j,:] = inf
            dists[:,k] = inf
            
        
    for i in arange(verticalized_parallels.shape[1]/2,verticalized_parallels.shape[1]-1,1):
        dists = zeros((verticalized_parallels.shape[0],len(parallels[i+1])))
        for j in range(verticalized_parallels.shape[0]):
            dists[j] = linalg.norm(array(parallels[i+1])[:,0:3]-verticalized_parallels[j][i][0:3],axis=1)
        for _ in range(len(dists)):
            j,k = unravel_index(argmin(dists),dists.shape)
            verticalized_parallels[j][i+1] = parallels[i+1][k]
            dists[j,:] = inf
            dists[:,k] = inf

    return verticalized_parallels
    
def compute_perpendicular_vector(vector_1):
    """
    Compute a perpendicular vector of a given vector.

    Keyword arguments:
    vector_1 -- given vector.
    """
    
    vector_1 = vector_1/linalg.norm(vector_1)
    vector_2 = copy(vector_1)
    while abs(dot(vector_1,vector_2)-1)<=1e-6:
        vector_2 = vector_1+random.uniform(1,-1,3)
        vector_2 = vector_2/linalg.norm(vector_2)                  
    return cross(vector_1,vector_2)

def distance_point_line_3d(x1, x2, x0):
    """
    A line in three dimensions specified by two points x1 and x2.
    The function calculates the euclidean distance between a point
    x0 and the line.
    
    keyword arguments:
    x1 -- (x,y,z) a point of the line
    x2 -- (x,y,z) a second point of the line
    x0 -- point not in the line to compute distance
    """

    return linalg.norm(cross(x0-x1,x0-x2))/linalg.norm(x2-x1)

def closest_point_line_3d(x1, x2, x0):
    """
    A line in three dimensions specified by two points x1 and x2.
    The function calculates the closest point in line to point
    x0.
    
    keyword arguments:
    x1 -- (x,y,z) a point of the line
    x2 -- (x,y,z) a second point of the line
    x0 -- point not in the line to compute distance
    """
    P10 = (x1-x0); P12 = (x1-x2)
    return -P12*dot(P10,P12)/dot(P12,P12)+x1

def isospherical_theta(xyz):
    return acos(xyz[2]/sqrt(xyz[0]**2+xyz[1]**2+xyz[2]**2))

def isospherical_phi(xyz):
    return atan2(xyz[1],xyz[0])

def isospherical_radius(xyz):
    return sqrt(xyz[0]**2+xyz[1]**2+xyz[2]**2)
    

def filter_by_distance(points, r = None, variation = 0.9, is_sorted = False):
    """
    The filter_by_distance method is an algorithm to delete the nearest neighbors
    points, inside a distance threshold. 

    Keyword arguments:
    points -- points to be filtered array(array).
    r -- distance threshold. If r is None, points are sorted w.r.t. x axis.
    variation -- min cossine of the angle between adjacent normals
    """

    if not is_sorted:
        points = points[argsort(points[:,0])]

    
    Tree = KDTree(points[:,0:3])
    rays = []
    N = len(points)
    I = ones((N, 1), dtype=bool)
    
    if r > 0 or r is not None:
        i=0
        while True:
            if I[i]:
                rays.append(points[i])
                idx = Tree.query_ball_point(points[i,0:3],r)
                idx = array(idx)
                idx = idx[idx>i]
                for j in idx:
                    if dot(points[i,3:6],points[j,3:6]) >= variation:
                        I[j]=False
            i+=1
            if i == N:
                break
    else:
        for i in range(N-1):
            if (dot(points[i,3:6],points[i-1,3:6]) < variation) or (dot(points[i,3:6],points[i+1,3:6]) < variation):
                rays += [points[i]]
        
    
    return array(rays)

class IterSurface:
    """ Inheritable class to surfaces that can be iterated and generate the coating
    trajectories.

    Keyword arguments:
    Rn0 -- initial parameter of the iteration
    stopR -- stop parameter of the iteration
    coatingstep -- iter step
    """

    __metaclass__ = ABCMeta
    
    def __init__(self, Rn0=3.770, stopR=1.59, coatingstep = 0.003):
        self._Rn0 = Rn0
        self._Rn = Rn0
        self.stopR = stopR
        self.coatingstep = coatingstep

    @abstractmethod
    def update(self):
        pass

    @abstractmethod
    def f(self, p):
        pass 

    @abstractmethod
    def df(self, p):
        pass

    @abstractmethod
    def f_array(self, p):
        pass 

    @abstractmethod
    def find_iter(self, p0):
        pass

    @abstractmethod
    def criteria(self):
        pass 

    @abstractmethod
    def findnextparallel(self, p):
        pass

    @abstractmethod
    def name(self, p):
        pass


class Sphere(IterSurface):
    """ Sphere surface class. An object of this class can be
    a surface to be iterated and generate the coating trajectories.

    Keyword arguments:
    Rn0 -- initial parameter of the sphere
    stopR -- stop parameter of the iteration
    coatingstep -- iter step
    """

    def __init__(self, Rn0=3.770, stopR=1.59, coatingstep = 0.003, center = [0,0,0]):
        IterSurface.__init__(self, Rn0, stopR, coatingstep)
        self.center = center

    def update(self):
        self._Rn = self._Rn-self.coatingstep
        
    def f(self, p):
        return ((p[0]-self.center[0])**2 +
                (p[1]-self.center[1])**2 +
                (p[2]-self.center[2])**2 -
                self._Rn**2)

    def df(self, p):
        return array([2*p[0]-self.center[0],
                      2*p[1]-self.center[1],
                      2*p[2]-self.center[2]])

    def f_array(self, p):
        if len(p)==1: return self.f(p[0])
        p = array(p)
        return sum((p[:,0:3]-self.center)*(p[:,0:3]-self.center),1)-self._Rn**2
        
    def find_iter(self, p0):
        self._Rn = sqrt((p0[0]-self.center[0])**2+
                        (p0[1]-self.center[1])**2+
                        (p0[2]-self.center[2])**2)
        return

    def criteria(self):
        return self._Rn>self.stopR

    def findnextparallel(self, p):
        d = sqrt((p[0]-self.center[0])**2+
                 (p[1]-self.center[1])**2+
                 (p[2]-self.center[2])**2)
        n = ceil((d-self.stopR)/self.coatingstep)
        self._Rn = min(self.stopR+n*self.coatingstep, self._Rn)
        return

    def name(self):
        return 'Sphere'

class Plane(IterSurface):
    """ Plane surface class. An object of this class can be
    a surface to be iterated and generate the coating trajectories.
    Keyword arguments:
    Rn0 -- height of the plane (z)
    stopR -- stop parameter of the iteration
    coatingstep -- iter step
    """

    def __init__(self, Rn0=1, stopR=-4, coatingstep = 0.003, k=0):
        IterSurface.__init__(self, Rn0, stopR, coatingstep)
        self.k = k

    def update(self):
        self._Rn = self._Rn-self.coatingstep

    def f(self, p):
        return p[self.k]-self._Rn

    def df(self, p):
        n = array([0, 0, 0])
        n[self.k] = 1
        return n

    def f_array(self, p):
        if len(p)==1: return self.f(p[0])
        p = array(p)
        return p[:,self.k]-self._Rn

    def find_iter(self, p0):
        self._Rn = p0[self.k]

    def criteria(self):
        return self._Rn>self.stopR

    def findnextparallel(self, p):
        d = p[self.k]
        n = ceil((d-self.stopR)/self.coatingstep)
        self._Rn = min(self.stopR+n*self.coatingstep, self._Rn)
        return

    def name(self):
        return 'Plane'


class Plane2(IterSurface):
    """ Plane surface class. An object of this class can be
    a surface to be iterated and generate the coating trajectories.
    plane: ax+by+cz+d=0
    normal_plane = [a,b,c]
    
    Keyword arguments:
    Rn0 -- d
    stopR -- stop parameter of the iteration
    coatingstep -- iter step
    normal_plane -- [a,b,c]
    """

    def __init__(self, Rn0=1, stopR=-4, coatingstep = 0.003, normal_plane=array([0,1,0])):
        IterSurface.__init__(self, Rn0, stopR, coatingstep)
        self.normal_plane = normal_plane

    def update(self):
        self._Rn = self._Rn-self.coatingstep

    def f(self, p):
        return (dot(self.normal_plane, p[0:3])+self._Rn)/linalg.norm(self.normal_plane)

    def df(self, p):
        return array(self.normal_plane)

    def f_array(self, p):
        if len(p)==1: return self.f(p[0])
        p = array(p)
        return (dot(p[:,0:3],self.normal_plane)+self._Rn)/linalg.norm(self.normal_plane)

    def find_iter(self, p0):
        self._Rn = -dot(self.normal_plane,p0[0:3])

    def criteria(self):
        return self._Rn>self.stopR

    def findnextparallel(self, p):
        self.find_iter(p)
        self.update()
        return

    def name(self):
        return 'Plane2'



_backdif2 = array([[1.0, -1.0],
                   [0,    0  ]])

_backdif3 = array([[3.0/2,    -2.0,   1.0/2 ],
                   [1.0,      -2.0,   1.0   ]])

_backdif4 = array([[11.0/6,   -3.0,   3.0/2,  -1.0/3],
                   [2.0,      -5.0,   4.0,    -1.0  ]])

_backdif5 = array([[25.0/12,  -4.0,    3.0,      -4.0/3,  1.0/4],
                   [35.0/12,  -26.0/3, 19.0/2,   -14.0/3, 11.0/12]])

_backdif6 = array([[137.0/60, -5.0,     5.0,      -10.0/3, 5.0/4,   -1.0/5],
                   [15.0/4,   -77.0/6,  107.0/6,  -13.0,   61.0/12, -5/6]])

_backdif7 = array([[49.0/20,  -6.0,     15.0/2,  -20.0/3,  15.0/4, -6.0/5,  1.0/6],
                   [203.0/45, -87.0/5,	117.0/4, -254.0/9, 33.0/2, -27.0/5, 137.0/180]])

_backdiflist = (_backdif2,
                _backdif3,
                _backdif4,
                _backdif5,
                _backdif6,
                _backdif7)

def backward_difference(turbine, joints_trajectory):
    """
    joints_trajectory - joints_trajectory[0] is the most recent
    """
    joints = array(joints_trajectory[:7])
    h = turbine.config.model.trajectory_step
    size = len(joints)

    w, alpha = h**(size-size/2) * dot(_backdiflist[size-2],
                                      h**(size/2) * joints)
    
    return w, alpha

def partial_backward_difference(turbine, old_joints_trajectory):
    """
    old_joints_trajectory - old_joints_trajectory[0] is the most recent past value
    Compute from old_joints_trajectory the partial backward diffeference, before adding the last point
    """
    joints = array(old_joints_trajectory[:6])
    h = turbine.config.model.trajectory_step
    size = len(joints)+1

    w_part, alpha_part = dot(_backdiflist[size-2][:,1:],
                                      h**(size/2) * joints)
    
    return w_part, alpha_part, size

def update_backward_difference(turbine, actual_joint, w_part, alpha_part, size):
    """
    joints_trajectory - joints_trajectory[0] is the most recent
    """
    h = turbine.config.model.trajectory_step

    w_alpha = h**(size-size/2)* (transpose((w_part, alpha_part)) +
                                  multiply.outer(h**(size/2) * actual_joint,
                                                 _backdiflist[size-2][:,0]))
        

    return w_alpha[:,:,0],w_alpha[:,:,1]

def general_finite_difference(time, joints, times):
    """
    Compute finite difference for uneven grid
    """
    
    if len(joints) < 3 or len(times) < 3:
        raise ValueError('Too short segment')
    if len(joints) != len(times):
        raise ValueError('Unequal length of times and joints')

    N = len(times)
    M=2
    s = zeros((M+1,N)) #the coeficients of finite difference for up to second derivative
    s_temp = array(s)
    s[0,0] = 1.
    
    dtN = 1. #maximum power of dt
    for n in range(1,N):
        dtn = 1. #powers of dt
        for v in range(n):
            dt = times[n]-times[v]
            dtn *= dt
            for m in range(min(n,M)+1):
                s_temp[m,v] = ((times[n]-time)*s[m,v] - m*s[m-1,v])/dt
                
        for m in range(min(n,M)+1):
            s_temp[m,n] = ( m*s[m-1,n-1] - (times[n-1]-time)*s[m,n-1])*(dtN/dtn)

        s = array(s_temp)

        dtN = dtn

    df = dot(s,joints)
    return df[0],df[1],df[2] #joint, w, alpha
    
def roundrobin(*iterables):
    "roundrobin('ABC', 'D', 'EF') --> A D E B F C"
    # Recipe credited to George Sakkis
    pending = len(iterables)
    nexts = cycle(iter(it).next for it in iterables)
    while pending:
        try:
            for next in nexts:
                yield next()
        except StopIteration:
            pending -= 1
            nexts = cycle(islice(nexts, pending))    

def intersection(points1, points2, r):
    intersections = []
    points1 = array(list(points1))
    points2 = array(list(points2))
    for point1 in points1:
        in_radius = points2[linalg.norm(points2-point1,axis=1)<=r]
        if len(in_radius)>0:
            intersections.extend(list((in_radius+point1)/2))
    intersections_tuple = [tuple(i) for i in intersections]
    return intersections_tuple

def distance_line_bases(x1, x2, bases, min_threshold, max_threshold):
    distance_str = 0; distance = None; point_near = []
    for base in bases:
        angle = sign(base[1])*atan2(x1[0]-base[0],abs(base[1]-x1[1]))
        if sign(x1[0])!=sign(base[0]):
            continue
        if angle<0 and abs(x1[0])<abs(base[0]):
            continue
        d = distance_point_line_3d(array(x1), array(x2), array(base))
        if d<=max_threshold:
            distance_str+=1
        if d<=min_threshold:
            if distance>d or distance==None:
                distance = d
                point_near = base
    return point_near, distance, distance_str
        
def union_line_grids(line_grid, lines):
    """
    Method return a set which is the union of given sets.

    Keyword arguments:
    line_grid -- it is a dictionary. line_grid relates lines (rays)
    with coatable grids for that line. This function makes sense
    with the DB pkl file.
    lines -- lines to unify grids.
    """
    
    line_union = set()
    for line in lines:
        line_union = line_union.union(line_grid[line])
    return line_union

def linear_interpolation_ray(ray0, ray1, threshold=1e-2):
    """
    Return the linear interpolation from ray0 to ray1
    """

    d = linalg.norm(ray0[0:3]-ray1[0:3])
    n = int(d/threshold)
    points = []
    for i in range(0,n):
        t = threshold*i/d
        p = (1-t)*ray0+ray1*t
        p[3:] = p[3:]/linalg.norm(p[3:])
        points.append(p)
    return points

def linear_interpolation_points(point0, point1, threshold=1e-2):
    """
    Return the linear interpolation from point0 to point1
    """

    d = linalg.norm(point0-point1)
    n = ceil(d*1./threshold)+1
    points = []
    lin = linspace(0,1,n)
    for t in lin:
        p = (1-t)*point0+point1*t
        points.append(p)
    return points, lin[1]*d

def linear_interpolation_joint(joint0, joint1, threshold=1e-4) :
    """
    Return the linear interpolation from joint0 to joint1
    """

    d = linalg.norm(joint0-joint1)
    n = ceil(d/threshold)+1
    joints = []
    for t in linspace(0,1,n):
        joint = (1-t)*joint0+joint1*t
        joints.append(joint)
    return joints
    
def homogenous_matrix_cubic_interpolation(Ai,Aj,Bi,Bj,n):
    t = linspace(0,1,n)
    T = []
    invAi = linalg.inv(Ai)
    invAj = linalg.inv(Aj)
    Ci = logm(einsum('ik,km,mj->ij',expm(Bi),invAi,Aj))
    Di = Bj - einsum('ik,kl,ln,nm,mj->ij',invAj,Ai,Bi,invAi,Aj) - 2*Ci
    for ti in t:
        T.append(dot(Ai,einsum('ik,km,mj->ij',
                               expm(ti*Bi), expm(ti**2*Ci), expm((ti**3-ti**2)*Di))))
    return T

def lenlist(list_of_lists):
    count = 0
    for l in list_of_lists:
        count+=len(l)
    return count

def notempty(l):
    return [x for x in l if x != []]

def makeTbyRP(R,P):
    T = eye(4,4)
    T[0:3,0:3]=R
    T[0:3,3]=P
    return T

def get_manip_transforms(robot, joints, rotation=False, position=False):
    T = []
    R = []
    p = []
    manip = robot.GetActiveManipulator()
    with robot:
        for joint in joints:
            robot.SetDOFValues(joint)
            Tee = manip.GetTransform()
            p.append(Tee[0:3,3])
            R.append(Tee[0:3,0:3])
            T.append(Tee)
    if rotation and position:
        return Tee, R, p
    if rotation:
        return Tee, R
    if position:
        return Tee, p
    return Tee
    
def filter_trajectory(points, threshold=1e-2, remove=False):
    new_points = copy(points)
    index_to_remove = []
    comp = points[0][0:3]
    for i in range(1,len(new_points)-1):
        d = linalg.norm(new_points[i][0:3]-comp)
        if d < threshold:
            index_to_remove.append(i)
        else: comp = new_points[i][0:3]

    if remove==True:
        to_remove = []
        for index in index_to_remove:
            to_remove.append(points[index])
    
    for index in sorted(index_to_remove, reverse=True):
        del new_points[index]

    if remove==True: return new_points, to_remove
    return new_points
    
