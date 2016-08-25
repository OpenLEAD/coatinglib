from numpy import load, savez_compressed, zeros, ones, arange, r_, c_, outer, tile
from numpy import meshgrid, array, shape, sum, eye, dot, argmin, concatenate, sqrt
from numpy import argmax
from os import makedirs
import errno
from openravepy import RaveCreateCollisionChecker, matrixFromAxisAngle
from openravepy.misc import SpaceSamplerExtra
from scipy.spatial import KDTree
import mathtools
from math import atan2, pi, ceil
from openrave_plotting import plot_points, plot_points_array, plot_point, remove_points
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

    def save_samples(self):
        savez_compressed('Blade/'+self._name+'_points.npz', array=self._points)
        return

    def save_model(self):
        savez_compressed('Blade/'+self._model.model_type+'/'+self._model._name+'_w.npz',
                        array=self._model._w)
        savez_compressed('Blade/'+self._model.model_type+'/'+self._model._name+'_points.npz',
                        array=self._model._points)
        return

    def save_trajectories(self):
        savez_compressed('Blade/Trajectory/'+self._model._name+'_trajectories.npz',
                                 array=self._trajectories)
        return

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

        if self.visualization:
            try:
                self.turbine.env.RemoveKinBody(self.turbine.primary)
            except: None
            try:
                self.turbine.env.RemoveKinBody(self.turbine.secondary)
            except: None
            try:
                self.turbine.env.RemoveKinBody(self.turbine.iris)
            except: None
            try:
                self.turbine.env.RemoveKinBody(self.turbine.runner_area)
            except: None
            try:
                self.turbine.env.RemoveKinBody(self.turbine.robot)
            except: None
            try:
                for i in range(1,len(self.turbine.blades)):
                    self.turbine.env.RemoveKinBody(self.turbine.blades[i])
            except: None
            
        
        Rminmax = [self.turbine.model.nose_radius-0.01, self.turbine.model.runner_radius+0.01]
     
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
        #self._points = self._points[sqrt(sum(self._points[:,0:3]*self._points[:,0:3],1))>Rminmax[0]]
        #self._points = self._points[sqrt(sum(self._points[:,0:3]*self._points[:,0:3],1))<Rminmax[1]]
        #self._points[:,0:3] = self._points[:,0:3] + 0.1*self._points[:,3:6] # Shifting points for tree filtering

        self._points = self.filter_by_distance(self._points, min_distance_between_points)
        #self._points[:,0:3] = self._points[:,0:3] - 0.1*self._points[:,3:6]
        self.save_samples()

        if self.visualization:
            plot_points(self.turbine, self._points, 'sampling', ((1,0,0)))

        return             
   
    def filter_by_distance(self, points, r):
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
        
    def make_model(self, iter_surface=None):
        """
        The make_model method is an algorithm to generate a mathematical representation
        of the object (mesh). This method can be called after the sampling method,
        and never before. The model is generated by the model_type (constructor), e.g. RBF.   

        Keyword arguments:
        iter_surface -- This is the surface to be iterated, as mathtools.sphere.
        If the model has more than 7000 points, this argument is needed.

        """
        print "BladeModeling::make_model - Warning: this is a data-intensive computing and might freeze your computer."

        if len(self._points)>7000:
            if iter_surface is None:
                raise ValueError("""The number of points is bigger than 7000. It is not safe
                                 to make an unique model that big. Create an iterate surface
                                 to partitionate the model (check mathtools.IterSurface).""")
            else:
                if not issubclass(iter_surface.__class__, mathtools.IterSurface):
                    raise IterSurfaceError()
                else:
                    points = self.split_blade_points(self._points, iter_surface)
                    
        else:
            self._model._points = self._points
            self._model.make()
            self.save_model()
        return


    def split_blade_points(self, points, iter_surface):
        #r0 = iter_surface._Rn
        #r_old = iter_surface._Rn
        #final_points = []
        #old_points = []
        number_of_parts = ceil(len(points)*1.0/4000)
        number_of_points_per_part = len(points)*1.0/number_of_parts
        #counter = 1
        points_distance = iter_surface.f_array(points)
        points = [x for (y,x) in sorted(zip(points_distance, points))]
        points = array(points)
        
        while iter_surface.criteria():
            r_new = iter_surface._Rn
            positive_values = iter_surface.f_array(points)>=0
            iter_surface._Rn = r_old
            negative_values = iter_surface.f_array(points)<=0
            iter_surface._Rn = r_new
            new_points = points[positive_values & negative_values]

            if counter!=number_of_parts:
                if len(new_points)<number_of_points_per_part:
                    old_points = deepcopy(new_points)
                else:
                    if len(old_points)==0:
                        raise ValueError('It is not possible to split the cloud with this surface/step.')
                    else:
                        print 'iter_surface._Rn = ', iter_surface._Rn
                        print 'r_old =', r_old
                        final_points.append(old_points)
                        old_points = []
                        r_old = iter_surface._Rn
                        counter += 1
            else:
                old_points = deepcopy(new_points)
                
            iter_surface.update()

        if len(old_points)!=0:
            final_points.append(old_points)

        iter_surface._Rn = r0    
        return final_points    
    
    

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
            P0 = copy(Pd)
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
            Rs = sqrt(sum(self._points[:,0:3]*self._points[:,0:3],1))
            Pd = self._points[argmax(Rs)]
            iter_surface.findnextparallel(Pd)
            Pd = mathtools.curvepoint(self._model, iter_surface, [Pd[0],Pd[1],Pd[2]])

        counter = 1
        while iter_surface.criteria():
            self._trajectories = drawParallel(self._trajectories, Pd, iter_surface)
            if counter%30==0:
                self.save_trajectories()
            p0=self._trajectories[-1][-1]
            iter_surface.update()
            Pd = mathtools.curvepoint(self._model, iter_surface, [p0[0],p0[1],p0[2]])
            counter+=1
        
        self.save_trajectories()
        return

    def remove_trajectories_from_env(self):
        remove_points(self.turbine, 'trajectories')

    def remove_sampling_from_env(self):
        remove_points(self.turbine, 'sampling')    
