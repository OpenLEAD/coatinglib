from numpy import load, savez_compressed, zeros, ones, arange, r_, c_, outer, tile
from numpy import meshgrid, array, shape, sum, eye, dot, argmin, concatenate, sqrt
from os import makedirs
import errno
from openravepy import RaveCreateCollisionChecker
from openravepy.misc import SpaceSamplerExtra
from scipy.spatial import KDTree
import mathtools
from math import atan2, pi
from openrave_plotting import plotPoints, plotPointsArray, plotPoint, removePoints
import copy

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
        self._modelLoaded = False
        self._samplingLoaded = False
        self.visualization = visualization
        if self.visualization:
            turbine.env.SetViewer('qtcoin')
        
        try:
            makedirs('./Blade')
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                print "Problem creating Blade folder"
                raise
        try:
            self._points = load('Blade/'+self._name+'_points.npz')
            self._points = self._points['array']
            print "BladeModeling::init - samples are loaded."
            self._samplingLoaded = True
        except:
            self._points = []
            print "BladeModeling::init - samples could not be loaded. Call method BladeModeling::sampling." 

        try:
            self._model._w = load('Blade/'+self._model.model_type+'/'+self._model._name+'_w.npz')
            self._model._w = self._model._w['array']
            self._model._points = load('Blade/'+self._model.model_type+'/'+self._model._name+'_points.npz')
            self._model._points = self._model._points['array']
            self._modelLoaded = True
            print "model::init - model is loaded."
        except:
            self._model._w = []
            self._model._points = self._points
            print "BladeModeling::init - model could not be loaded."

        try:
            makedirs('./Blade/Trajectory')
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                print "BladeModeling::init - Problem creating Blade/Trajectory folder"
                raise    
        try:
            self._trajectories = load('Blade/Trajectory/'+self._name+'_trajectories.npz')
            self._trajectories = self._trajectories['array']
            self._trajLoaded = True
            print "BladeModeling::init - Trajectories are loaded."
        except:
            self._trajectories = []
            self._trajLoaded = False
            print "BladeModeling::init - Trajectories could not be loaded."
            

    def sampling(self, delta = 0.005):
        print 'Blade::sampling - Warning: this is a data-intensive computing and might freeze your computer.'
        
        Rminmax = [self.turbine.model.nose_radius-0.01, self.turbine.model.runner_radius+0.01]
        bladerotation=[self.turbine.environment.blade_angle, 'y']
        while True:
            if self._samplingLoaded:
                answer = raw_input("You are performing resampling, as the samples were loaded. Are you sure you want to continue? [y or n]")
                if answer == 'y':
                    self._points = []
                    break
                elif answer=='n':return
            else: break    
            
        cc = RaveCreateCollisionChecker(self.turbine.env,'ode')
        if cc is not None:
                ccold = self.turbine.env.GetCollisionChecker()
                self.turbine.env.SetCollisionChecker(cc)
                cc = ccold
        self._blade.SetTransform(eye(4))
        self._blade = mathtools.Rotate(self._blade, -bladerotation[0], bladerotation[1])
        ab = self._blade.ComputeAABB()
        p = ab.pos()
        e = ab.extents()+0.01 # increase since origin of ray should be outside of object
        sides = array((
                     (e[0],0,0,-1,0,0,0,e[1],0,0,0,e[2]), #x
                     (-e[0],0,0,1,0,0,0,e[1],0,0,0,e[2]), #-x
                     (0,0,e[2],0,0,-1,e[0],0,0,0,e[1],0), #z
                     (0,0,-e[2],0,0,1,e[0],0,0,0,e[1],0), #-z
                     (0,e[1],0,0,-1,0,e[0],0,0,0,0,e[2])  #y
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
        def treeFilter(points, r=0.05):
            print "Blade::_treeFilter - Starting filtering points"
            T = KDTree(points[:,0:3])
            rays = []
            N = len(points)
            i=0
            I = ones((len(points), 1), dtype=bool)
            while True:
                if I[i]:
                    rays.append(points[i])
                    idx = T.query_ball_point(points[i,0:3],r)
                    idx = array(idx)
                    idx = idx[idx>i]
                    for j in idx:I[j]=False
                i+=1
                if i==len(points):break
            return array(rays)
        self._points = treeFilter(self._points)
        #self._points[:,0:3] = self._points[:,0:3] - 0.1*self._points[:,3:6]
        savez_compressed('Blade/'+self._name+'_points.npz', array=self._points)
        print "BladeModeling::samplig - terminates."
        self._samplingLoaded = True
        
        if self.visualization:
            self.turbine.env.RemoveKinBody(self.turbine.primary)  
            self.turbine.env.RemoveKinBody(self.turbine.secondary)
            plotPoints(self.turbine, self._points, 'sampling', ((1,0,0)))   
        
    def make_model(self):
        try:
            makedirs('./Blade/'+self._model.model_type)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                print "BladeModeling::make_model - Problem creating Blade/mode_type folder"
                raise
            
        self._model._points = self._points
        self._model.make()
        savez_compressed('Blade/'+self._model.model_type+'/'+self._model._name+'_w.npz', array=self._model._w)
        savez_compressed('Blade/'+self._model.model_type+'/'+self._model._name+'_points.npz', array=self._model._points)
        self._modelLoaded = True
        print "BladeModeling::make_model - terminates."
            

    def generate_trajectory(self, iter_surface):
        """ Method generate the coating trajectories. The trajectories are
        the intersection between two surfaces: the blade model, and the surface
        to be iterated, e.g. spheres. The algorithm
        follows the marching method, documentation available in:
        http://www.mathematik.tu-darmstadt.de/~ehartmann/cdgen0104.pdf
        page 94 intersection surface - surface.
        """
        print "BladeModeling::generate_trajectory - Warning: this is a data-intensive computing and might freeze your computer."

        if self.visualization:
            try:
                self.turbine.env.RemoveKinBody(self.turbine.primary)  
                self.turbine.env.RemoveKinBody(self.turbine.secondary)
            except: None
            
        self._blade.SetTransform(eye(4))
        self._blade = mathtools.Rotate(self._blade, -self.turbine.environment.blade_angle, 'y')
        
        if self._modelLoaded: None
        else:
            print "BladeModeling::generate_trajectory - Model is not loaded. Load the model first with make_model method"
            return

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
                    plotPoint(self.turbine, P, 'trajectories', ((0,0,1)))

                y.append(P)
        if self._trajLoaded:
            tempY = []
            for y in self._trajectories:
                tempY.append(list(y))   
            self._trajectories = tempY
            Pd = self._trajectories[-1][-1]
            Rn = iter_surface.find_iter(Pd)
            Pd = mathtools.curvepoint(self._model, iter_surface, [Pd[0],Pd[1],Pd[2]])
            if self.visualization:
                plotPointsArray(self.turbine, self._trajectories, 'trajectories', ((0,0,1)))
        else:
            Pd = self._points[argmin(self._points[:,1])]
            iter_surface.findnextparallel(Pd)
            Pd = mathtools.curvepoint(self._model, iter_surface, [Pd[0],Pd[1],Pd[2]])

        counter = 0
        while iter_surface.criteria():
            self._trajectories = drawParallel(self._trajectories, Pd, iter_surface)
            if counter%50==0:
                savez_compressed('Blade/Trajectory/'+self._name+'_trajectories.npz', array=self._trajectories)   
            p0=self._trajectories[-1][-1]
            iter_surface.update()
            Pd = mathtools.curvepoint(self._model, iter_surface, [p0[0],p0[1],p0[2]])
            counter+=1
        print "BladeModeling::generate_trajectory - terminates."
        savez_compressed('Blade/Trajectory/'+self._name+'_trajectories.npz', array=self._trajectories)

    def RemoveTrajectoriesFromEnv(self):
        removePoints(self.turbine, 'trajectories')

    def RemoveSamplingFromEnv(self):
        removePoints(self.turbine, 'sampling')    
