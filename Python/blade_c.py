from numpy import load, savez_compressed, zeros, ones, arange, r_, c_, outer, tile
from numpy import meshgrid, array, shape, sum, eye, dot, argmin, concatenate, sqrt
from os import makedirs
import errno
from openravepy import RaveCreateCollisionChecker
from openravepy.misc import SpaceSamplerExtra
from scipy.spatial import KDTree
import mathtools
from math import atan2, pi

class Blade:
    """ Blade class.

    Keyword arguments:
    name -- the name of the blade.
    env -- environment object.
    blade_string -- name of the blade in xml.
    """

    def __init__(self, name, env, blade_string):
        blade = next(body for body in env.GetBodies() if body.GetName()==blade_string)
        self._name = name
        self.blade = blade


class BladeModeling:
    """ BladeModeling class for blade modelling.

    Keyword arguments:
    name -- the name of the blade.
    model_type -- model object (e.g. RBF).
    env -- environment object
    blade -- body blade object.
    """

    def __init__(self, name, model_type, env, blade):
        self.handles = []
        self._name = name
        self._model = model_type
        self._env = env
        self._blade = blade
        self._modelLoaded = False
        
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

            

    def sampling(self, Rminmax=[1.59,3.75], delta = 0.01, coatingdistance = 0.23, bladerotation=[0,'']):
        print 'Balde::sampling - Warning: this is a data-intensive computing and might freeze your computer.'
        cc = RaveCreateCollisionChecker(self._env,'ode')
        if cc is not None:
                ccold = self._env.GetCollisionChecker()
                self._env.SetCollisionChecker(cc)
                cc = ccold
        self._blade.SetTransform(eye(4))
        self._blade = mathtools.Rotate(self._blade, bladerotation[0], bladerotation[1])
        ab = self._blade.ComputeAABB()
        p = ab.pos()
        e = ab.extents()+0.01 # increase since origin of ray should be outside of object
        sides = array((
                     (e[0],0,0,-1,0,0,0,e[1],0,0,0,e[2]),
                     (-e[0],0,0,1,0,0,0,e[1],0,0,0,e[2]),
                     (0,0,e[2],0,0,-1,e[0],0,0,0,e[1],0),
                     (0,0,-e[2],0,0,1,e[0],0,0,0,e[1],0)
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
            collision, info = self._env.CheckCollisionRays(rays,self._blade)
            # make sure all normals are the correct sign: pointing outward from the object)
            newinfo = info[collision,:]
            if len(newinfo) > 0:
                  newinfo[sum(rays[collision,3:6]*newinfo[:,3:6],1)>0,3:6] *= -1
                  self._points = r_[self._points,newinfo]
        self._points = self._points[sqrt(sum(self._points[:,0:3]*self._points[:,0:3],1))>Rminmax[0]]
        self._points = self._points[sqrt(sum(self._points[:,0:3]*self._points[:,0:3],1))<Rminmax[1]]
        self._points[:,0:3] = self._points[:,0:3] + coatingdistance*self._points[:,3:6]
        T = mathtools.T(bladerotation[0],bladerotation[1])
        self._points[:,0:3] = dot(self._points[:,0:3],T[0:3,0:3])
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
        savez_compressed('Blade/'+self._name+'_points.npz', array=self._points)
        print "BladeModeling::samplig - terminates."
        
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
        to be iterated, e.g. spheres (radius 1.425 to 3.770). The algorithm
        follows the marching method, documentation available in:
        http://www.mathematik.tu-darmstadt.de/~ehartmann/cdgen0104.pdf
        page 94 intersection surface - surface.
        """
        print "BladeModeling::generate_trajectory - Warning: this is a data-intensive computing and might freeze your computer."

        if self._modelLoaded: None
        else:
            print "BladeModeling::generate_trajectory - Model is not loaded. Load the model first with make_model method"
            return

        def plotPoint(env, point, handles, color):
            handles.append(env.plot3(points=array(point)[0:3], pointsize=5, colors=color))
            return handles

        def drawParallel(Y, Pd, iter_surface):
            dt = 3e-3
            theta0 = 180*atan2(-Pd[2],Pd[0])/pi
            counter = 0
            y = [Pd]       
            while True:
                tan = mathtools.surfaces_tangent(y[-1], iter_surface)
                P = mathtools.curvepoint(self._model, iter_surface, y[-1][0:3]-tan*dt)
                theta = 180*atan2(-P[2],P[0])/pi
                if counter==0:
                    if theta>theta0:
                        counter+=1
                else:
                    if theta<theta0:
                        Y.append(y)
                        return Y
                y.append(P)
                self.handles=plotPoint(self._env, P, self.handles, array((1,0,0)))

        try:
            makedirs('./Blade/Trajectory')
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                print "BladeModeling::generate_trajectory - Problem creating Blade/Trajectory folder"
                raise
        try:
            Y = load('Blade/Trajectory/'+self._name+'_trajectories.npz')
            Y = Y['array']
            tempY = []
            for y in Y:
                tempY.append(list(y))   
            Y = tempY
            Pd = Y[-1][-1]
            Rn = iter_surface.find_iter(Pd)
            Pd = mathtools.curvepoint(self._model, iter_surface, [Pd[0],Pd[1],Pd[2]])
            print "BladeModeling::generate_trajectory - Trajectories are loaded."
        except:
            print "BladeModeling::generate_trajectory - Trajectories could not be loaded."
            Y=[[]]
            Pd = self._points[argmin(self._points[:,1])]
            iter_surface.findnextparallel(Pd)
            Pd = mathtools.curvepoint(self._model, iter_surface, [Pd[0],Pd[1],Pd[2]])

        while iter_surface.criteria:
            Y = drawParallel(Y, Pd, iter_surface)
            savez_compressed('Blade/trajectory/'+self._name+'_trajectories.npz', array=Y)   
            p0=Y[-1][-1]
            iter_surface.update()
            Pd = mathtools.curvepoint(self._model, iter_surface, [p0[0],p0[1],p0[2]])    
        print "BladeModeling::generate_trajectory - terminates."
