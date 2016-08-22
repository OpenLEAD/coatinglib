from numpy import load, savez_compressed, zeros, ones, arange, r_, c_, outer, tile
from numpy import meshgrid, array, shape, sum, eye, dot, argmin, concatenate, sqrt
from os import makedirs
import errno
from openravepy import RaveCreateCollisionChecker, matrixFromAxisAngle
from openravepy.misc import SpaceSamplerExtra
from scipy.spatial import KDTree
import mathtools
from math import atan2, pi
from openrave_plotting import plot_points, plot_points_array, plot_point, remove_points
import copy

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

    def __init__(self, name, model_type, turbine, visualization = False):
        self.turbine = turbine
        self._name = name
        self._model = model_type
        self._blade = turbine.blades[0]
        self._trajectories = []
        self._points = []
        self._model._w = []
        self._trajectories = []
        self._model._points = self._points
        self.visualization = visualization

        if self.visualization:
            turbine.env.SetViewer('qtcoin')
        
        try:
            makedirs('./Blade')
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
        try:
            makedirs('./Blade/'+self._model.model_type)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
        try:
            makedirs('./Blade/Trajectory')
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise    

    def load_samples(self):
        try:
            self._points = load('Blade/'+self._name+'_points.npz')
            self._points = self._points['array']
            print "Samples are loaded."
            self._samplingLoaded = True
        except IOError:
            raise('Samples could not be loaded. Call method BladeModeling::sampling')
        return

    def load_model(self):
        try:
            self._model._w = load('Blade/'+self._model.model_type+'/'+self._model._name+'_w.npz')
            self._model._w = self._model._w['array']
            self._model._points = load('Blade/'+self._model.model_type+'/'+self._model._name+'_points.npz')
            self._model._points = self._model._points['array']
            self._modelLoaded = True
            print "Model is loaded."
        except IOError:
            raise("Model could not be loaded.")
        return
    
    def load_trajectories(self):
        try:
            self._trajectories = load('Blade/Trajectory/'+self._model._name+'_trajectories.npz')
            self._trajectories = self._trajectories['array']
            print "Trajectories are loaded."
        except IOError:
            raise("Trajectories could not be loaded.")

    def sampling(self, delta = 0.005, min_distance_between_points=0.05):
        """ The sampling method is an algorithm for uniformly sampling objects.
        It computes the bounding box of the object (object inscribed by cube), uniformly
        samples the box's faces, and checks rays collision from cube's samples to the object.
        Rays with collision and its normal vectors are stored as object's samples.
        Finally, a kdtree is computed from object's samples and neighbooring points are removed  
        to generate an uniformly sampling on the object (filtering data by distance). 

        Keyword arguments:
        delta -- uniform distance between the cube's samples.
        min_distance_between_points -- uniform distance between the object's samples.
        min_distance_between_points >= delta
        """
        
        print 'Blade::sampling - Warning: this is a data-intensive computing and might freeze your computer.'
        
        Rminmax = [self.turbine.model.nose_radius-0.01, self.turbine.model.runner_radius+0.01]
        bladerotation=[self.turbine.environment.blade_angle, 'z']
     
        cc = RaveCreateCollisionChecker(self.turbine.env,'ode')
        if cc is not None:
                ccold = self.turbine.env.GetCollisionChecker()
                self.turbine.env.SetCollisionChecker(cc)
                cc = ccold
        self._blade.SetTransform(eye(4))
        self._blade.SetTransform(dot(self._blade.GetTransform(),
                                     matrixFromAxisAngle([0, -self.turbine.environment.blade_angle, 0]))
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
        self._points = self._points[sqrt(sum(self._points[:,0:3]*self._points[:,0:3],1))>Rminmax[0]]
        self._points = self._points[sqrt(sum(self._points[:,0:3]*self._points[:,0:3],1))<Rminmax[1]]
        #self._points = self._points[dot(self._points[:,3:6],[0,1,0])<0.8] # Filtering normals close to [0,1,0]
        #self._points[:,0:3] = self._points[:,0:3] + 0.1*self._points[:,3:6] # Shifting points for tree filtering
        def treeFilter(points, r):
            print "Blade::_treeFilter - Starting filtering points"
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
        self._points = treeFilter(self._points, min_distance_between_points)
        #self._points[:,0:3] = self._points[:,0:3] - 0.1*self._points[:,3:6]
        savez_compressed('Blade/'+self._name+'_points.npz', array=self._points)

        if self.visualization:
            self.turbine.env.RemoveKinBody(self.turbine.primary)  
            self.turbine.env.RemoveKinBody(self.turbine.secondary)
            plot_points(self.turbine, self._points, 'sampling', ((1,0,0)))
   
        
    def make_model(self):
        """
        The make_model method is an algorithm to generate a mathematical representation
        of the object (mesh). This method can be called after the sampling method,
        and never before. The model is generated by the model_type (constructor), e.g. RBF.   
        """
        print "BladeModeling::make_model - Warning: this is a data-intensive computing and might freeze your computer."
        self._model._points = self._points
        self._model.make()
        savez_compressed('Blade/'+self._model.model_type+'/'+self._model._name+'_w.npz',
                        array=self._model._w)
        savez_compressed('Blade/'+self._model.model_type+'/'+self._model._name+'_points.npz',
                        array=self._model._points)       

    def generate_trajectory(self, iter_surface):
        """
        Method generate the coating trajectories. The trajectories are
        the intersection between two surfaces: the blade model, and the surface
        to be iterated, e.g. spheres. The algorithm
        follows the marching method, documentation available in:
        http://www.mathematik.tu-darmstadt.de/~ehartmann/cdgen0104.pdf
        page 94 intersection surface - surface.

        Keyword arguments:
        iter_surface -- surface to be iterated, as mathtools.sphere.
        """
        print "BladeModeling::generate_trajectory - Warning: this is a data-intensive computing and might freeze your computer."

        if self.visualization:
            try:
                self.turbine.env.RemoveKinBody(self.turbine.primary)  
                self.turbine.env.RemoveKinBody(self.turbine.secondary)
            except: None
            
        self._blade.SetTransform(eye(4))
        self._blade.SetTransform(dot(self._blade.GetTransform(),
                                     matrixFromAxisAngle([0, -self.turbine.environment.blade_angle, 0]))
                                 )
        if len(self._model._w)>0: None
        else:
            print "BladeModeling::generate_trajectory - Model is not loaded. Load the model first with make_model method"
            return

        if not issubclass(iter_surface.__class__, mathtools.IterSurface):
            raise IterSurfaceError()

        def drawParallel(Y, Pd, iter_surface):
            dt = 3e-3
            P0 = copy.copy(Pd)
            counter = 0
            y = [Pd]
            while True:
                tan = mathtools.surfaces_tangent(y[-1], iter_surface)
                P = mathtools.curvepoint(self._model, iter_surface, y[-1][0:3]-tan*dt)
                if counter==0:
                    dP = abs(P-P0)
                    if dP[0]<=dt:
                        if dP[1]<=dt:
                            if dP[2]<=dt:
                                counter+=1
                else:
                    dP = abs(P-P0)
                    if dP[0]<=dt:
                        if dP[1]<=dt:
                            if dP[2]<=dt:
                                Y.append(y)
                                return Y
                if self.visualization:
                    plot_point(self.turbine, P, 'trajectories', ((0,0,0)))
                y.append(P)
                
        if len(self._trajectories)>0:
            tempY = []
            for y in self._trajectories:
                tempY.append(list(y))   
            self._trajectories = tempY
            Pd = self._trajectories[-1][-1]
            Rn = iter_surface.find_iter(Pd)
            Pd = mathtools.curvepoint(self._model, iter_surface, [Pd[0],Pd[1],Pd[2]])
            if self.visualization:
                plot_points_array(self.turbine, self._trajectories,
                                  'trajectories', ((0,0,0)))
        else:
            Pd = self._points[argmin(self._points[:,2])]
            iter_surface.findnextparallel(Pd)
            Pd = mathtools.curvepoint(self._model, iter_surface, [Pd[0],Pd[1],Pd[2]])

        counter = 1
        while iter_surface.criteria():
            self._trajectories = drawParallel(self._trajectories, Pd, iter_surface)
            if counter%30==0:
                savez_compressed('Blade/Trajectory/'+self._model._name+'_trajectories.npz',
                                 array=self._trajectories)
            p0=self._trajectories[-1][-1]
            iter_surface.update()
            Pd = mathtools.curvepoint(self._model, iter_surface, [p0[0],p0[1],p0[2]])
            counter+=1
        
        savez_compressed('Blade/Trajectory/'+self._model._name+'_trajectories.npz',
                         array=self._trajectories)

    def remove_trajectories_from_env(self):
        remove_points(self.turbine, 'trajectories')

    def remove_sampling_from_env(self):
        remove_points(self.turbine, 'sampling')    
