from numpy import load, savez_compressed, zeros, ones, arange, r_, c_, outer, tile
from numpy import meshgrid, array, shape, sum, eye, dot, argsort, concatenate, sqrt
from numpy import argmax, argmin
from os import makedirs
import errno
from openravepy import RaveCreateCollisionChecker, matrixFromAxisAngle
from scipy.spatial import KDTree
import mathtools
from math import pi, ceil
from copy import copy, deepcopy

class IterSurfaceError(Exception):    
    def __init__(self):
        Exception.__init__(self, "object is not a valid surface.")

class BladeModeling:
    """ BladeModeling class for blade modelling.

    Keyword arguments:
    name -- the name of the blade.
    model_type -- model object (e.g. RBF).
    env -- environment object
    blade -- body blade object.
    """

    turbine = None
    name = None
    blade = None
    _trajectories = []
    _points = []
    
    def __init__(self, name, turbine, blade):
        self.turbine = turbine
        self._name = name
        self._blade = blade    

    def save_samples(self):
        try:
            makedirs('./Blade')
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise   
        savez_compressed('Blade/'+self._name+'_points.npz', array=self._points)
        return

    def save_model(self, model):
        try:
            makedirs('./Blade/'+model.model_type)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise  
        savez_compressed('Blade/'+model.model_type+'/'+model._name+'_w.npz',
                        array=model._w)
        savez_compressed('Blade/'+model.model_type+'/'+model._name+'_points.npz',
                        array=model._points)
        return

    def save_trajectories(self, model):
        try:
            makedirs('./Blade/Trajectory')
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
        savez_compressed('Blade/Trajectory/'+model._name+'_trajectories.npz',
                                 array=self._trajectories)
        return

    def load_samples(self):
        try:
            self._points = load('Blade/'+self._name+'_points.npz')
            self._points = self._points['array']
        except IOError:
            raise IOError('Samples could not be loaded. Call method BladeModeling::sampling')
        return

    def load_model(self, model):
        try:
            model._w = load('Blade/'+model.model_type+'/'+model._name+'_w.npz')
            model._w = model._w['array']
            model._points = load('Blade/'+model.model_type+'/'+model._name+'_points.npz')
            model._points = model._points['array']
        except IOError:
            raise IOError("Model could not be loaded.")
        return
    
    def load_trajectories(self, model):
        try:
            self._trajectories = load('Blade/Trajectory/'+model._name+'_trajectories.npz')
            self._trajectories = self._trajectories['array']
            self._trajectories = self._trajectories.tolist()
        except IOError:
            raise IOError("Trajectories could not be loaded.")
        return

    def sampling(self, delta = 0.005, min_distance_between_points=0.05):
        """ The sampling method is an algorithm for uniformly sampling objects.
        It computes the bounding box of the object (object inscribed by cube), uniformly
        samples the box's faces, and checks rays collision from cube's samples to the object.
        Rays with collision and its normal vectors are stored as object's samples.
        Finally, a kdtree is computed from object's samples and neighbooring points are removed  
        to generate an uniformly sampling on the object (filtering data by distance).

        This is a data-intensive computing and might freeze your computer.

        Keyword arguments:
        delta -- uniform distance between the cube's samples.
        min_distance_between_points -- uniform distance between the object's samples.
        min_distance_between_points >= delta
        """
        
        cc = RaveCreateCollisionChecker(self.turbine.env,'ode')
        if cc is not None:
                ccold = self.turbine.env.GetCollisionChecker()
                self.turbine.env.SetCollisionChecker(cc)
                cc = ccold
 
        self._blade.SetTransform(eye(4))
        self._blade.SetTransform(dot(self._blade.GetTransform(),
                                     matrixFromAxisAngle([0, -self.turbine.config.environment.blade_angle, 0]))
                                 )
        
        ab = self._blade.ComputeAABB()
        p = ab.pos()
        e = ab.extents()+0.01 # increase since origin of ray should be outside of object
        sides = array((
                     (e[0],0,0,-1,0,0,0,e[1],0,0,0,e[2]), #x
                     (-e[0],0,0,1,0,0,0,e[1],0,0,0,e[2]), #-x
                     (0,0,e[2],0,0,-1,e[0],0,0,0,e[1],0), #z
                     (0,0,-e[2],0,0,1,e[0],0,0,0,e[1],0), #-z
                     (0,e[1],0,0,-1,0,e[0],0,0,0,0,e[2]), #y
                     (0,-e[1],0,0,1,0,e[0],0,0,0,0,e[2])  #-y
                     ))
        maxlen = 2*sqrt(sum(e**2))+0.03
        self._points = zeros((0,6))
        for side in sides:
            ex = sqrt(sum(side[6:9]**2))
            ey = sqrt(sum(side[9:12]**2))
            XX,YY = meshgrid(r_[arange(-ex,-0.25*delta,delta),0,arange(delta,ex,delta)],
                                  r_[arange(-ey,-0.25*delta,delta),0,arange(delta,ey,delta)])
            localpos = outer(XX.flatten(),side[6:9]/ex)+outer(YY.flatten(),side[9:12]/ey)
            N = localpos.shape[0]
            rays = c_[tile(p+side[0:3],(N,1))+localpos,maxlen*tile(side[3:6],(N,1))]
            collision, info = self.turbine.env.CheckCollisionRays(rays, self._blade)
            # make sure all normals are the correct sign: pointing outward from the object)
            newinfo = info[collision,:]
            if len(newinfo) > 0:
                  newinfo[sum(rays[collision,3:6]*newinfo[:,3:6],1)>0,3:6] *= -1
                  self._points = r_[self._points,newinfo]

        self._points = self.filter_by_distance(self._points, min_distance_between_points)
        self.save_samples()
        return             
   
    def filter_by_distance(self, points, r):
        """
        The filter_by_distance method is an algorithm to delete the nearest neighbors
        points, inside a distance threshold. 

        Keyword arguments:
        points -- points to be filtered array(array).
        r -- distance threshold. If r is None, points are sorted w.r.t. x axis.
        """

        points = points[argsort(points[:,0])]
        if r is None:
            return points
        
        Tree = KDTree(points[:,0:3])
        rays = []
        N = len(points)
        i=0
        I = ones((len(points), 1), dtype=bool)
        while True:
            if I[i]:
                rays.append(points[i])
                idx = Tree.query_ball_point(points[i,0:3],r)
                idx = array(idx)
                idx = idx[idx>i]
                for j in idx:I[j]=False
            i+=1
            if i==len(points):break
        return array(rays)
        
    def make_model(self, model, iter_surface=None):
        """
        The make_model method is an algorithm to generate a mathematical representation
        of the object (mesh). This method can be called after the sampling method,
        and never before. The model is generated by the model_type (constructor), e.g. RBF.

        This is a data-intensive computing and might freeze your computer.

        Keyword arguments:
        iter_surface -- This is the surface to be iterated, as mathtools.sphere.
        If the model has more than 9000 points, this argument is needed.
        """
        
        if len(self._points)>9000:
            if iter_surface is None:
                raise ValueError("The number of points is "+str(len(self._points))+""", which is
                                bigger than 9000. It is not safe to make an unique model this big.
                                Create an iterate surface to partitionate the model
                                (check mathtools.IterSurface).""")
            elif not issubclass(iter_surface.__class__, mathtools.IterSurface):
                    raise IterSurfaceError()
            else:
                raise ValueError("Data is too big")
                #points = self.split_blade_points(self._points, iter_surface)
                    
        else:
            model._points = self._points
            model.make()
            self.save_model(model)
        return

    def compute_initial_point(self, model, iter_surface):
        """
        The compute_initial_point computes the initial point to start the generating trajectories
        algorithm. If trajectories were loaded, the initial point is the last computed point projected
        Keyword arguments:
        point_on_surfaces -- initial point on both surfaces.
        iter_surface -- surface to be iterated, as mathtools.sphere.
        step -- it must be small, e.g. 1e-3. Otherwise the method will fail.
        """
                
        if len(self._trajectories)>0:
            last_computed_point = self._trajectories[-1][-1]
            iter_surface.find_iter(last_computed_point)
            iter_surface.update()
            point_on_surfaces = mathtools.curvepoint(model, iter_surface, last_computed_point[0:3])
        else:    
            initial_point = self._points[argmax(iter_surface.f_array(self._points))]
            iter_surface.findnextparallel(initial_point)
            point_on_surfaces = mathtools.curvepoint(model, iter_surface, initial_point[0:3])
        return point_on_surfaces

    def draw_parallel(self, point_on_surfaces, model, iter_surface, step):
        """
        The draw_parallel generates one coating trajectory. The trajectory is
        the intersection between two surfaces: the blade model, and the surface
        to be iterated, e.g. spheres. The algorithm
        follows the marching method, documentation available in:
        http://www.mathematik.tu-darmstadt.de/~ehartmann/cdgen0104.pdf
        page 94 intersection surface - surface.
        The function stops when a closed curve is found, i.e., when it is found a
        point near the initial point (<step distance).

        This is a data-intensive computing and might freeze your computer.

        Keyword arguments:
        point_on_surfaces -- initial point on both surfaces.
        iter_surface -- surface to be iterated, as mathtools.sphere.
        step -- it must be small, e.g. 1e-3. Otherwise the method will fail.
        """
        
        initial_point = copy(point_on_surfaces)
        counter = 0
        trajectory = [point_on_surfaces]
        while True:
            tan = mathtools.surfaces_tangent(trajectory[-1], iter_surface)
            next_point_on_surfaces = mathtools.curvepoint(model,
                                                          iter_surface,
                                                          trajectory[-1][0:3]-tan*step)
            dP = abs(next_point_on_surfaces[0:3]-initial_point[0:3])
            if counter==0:
                if max(dP)<=step:counter+=1
            else:
                if max(dP)<=step:
                    try:
                        p0 = trajectory[0][0:3]; p1 = trajectory[1][0:3]; p2 = trajectory[2][0:3]
                        if (max(abs(p0-p1))<=step) and (max(abs(p0-p2))>=step):
                            trajectory.append(next_point_on_surfaces)
                            return trajectory
                    except IndexError:
                        raise IndexError('Step is too big and the function terminated soon.')
            trajectory.append(next_point_on_surfaces)
        return   

    def generate_trajectories(self, model, iter_surface, step = 1e-3):
        """
        Method generate the coating trajectories. The trajectories are
        the intersection between two surfaces: the blade model, and the surface
        to be iterated, e.g. spheres. The algorithm
        follows the marching method, documentation available in:
        http://www.mathematik.tu-darmstadt.de/~ehartmann/cdgen0104.pdf
        page 94 intersection surface - surface.

        This is a data-intensive computing and might freeze your computer.

        Keyword arguments:
        iter_surface -- surface to be iterated, as mathtools.sphere.
        step -- it must be small, e.g. 1e-3. Otherwise the method will fail.
        """

        if not issubclass(iter_surface.__class__, mathtools.IterSurface):
            raise IterSurfaceError()

        self._blade.SetTransform(eye(4))
        self._blade.SetTransform(dot(self._blade.GetTransform(),
                                     matrixFromAxisAngle([0, -self.turbine.config.environment.blade_angle, 0]))
                                 )
        point_on_surfaces = self.compute_initial_point(model, iter_surface)

        try: 
            while iter_surface.criteria():
                self._trajectories.append(self.draw_parallel(point_on_surfaces, model, iter_surface, step))
                p0=self._trajectories[-1][-1]
                iter_surface.update()
                point_on_surfaces = mathtools.curvepoint(model, iter_surface, p0[0:3])
        except KeyboardInterrupt:  
            self.save_trajectories(model)
            raise
        self.save_trajectories(model)
        return 
