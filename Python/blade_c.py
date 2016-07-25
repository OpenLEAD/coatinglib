from numpy import load, savez_compressed, zeros, ones, arange, r_, c_, outer, tile, meshgrid, array, shape, sum, eye, dot
from os import makedirs
import errno
from openravepy import RaveCreateCollisionChecker
from openravepy.misc import SpaceSamplerExtra
from scipy.spatial import KDTree
import mathtools
from math import sqrt

class Blade:
    """ Blade class.

    Keyword arguments:
    name -- the name of the blade.
    env -- environment object
    """

    def __init__(self, name, env):
        blade = next(body for body in env.GetBodies() if body.GetName()=='pa1')
        self._name = name
        self.blade = blade


class BladeModeling:
    """ BladeModeling class for blade modelling.

    Keyword arguments:
    name -- the name of the blade.
    model_type -- RBF object.
    env -- environment object
    blade -- body blade object.
    """

    def __init__(self, name, model_type, env, blade):
        self._name = name
        self._model = model_type
        self._env = env
        bodies = env.GetBodies()
        self._blade = blade
        try:
            makedirs('./Blade')
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                print "Problem creating Blade folder"
                raise
        try:
            self._points = load('Blade/'+self._name+'_points.npz')
            self._points = self._points['array']
            print "A blade sampling was loaded."
        except:
            self._points = []
            print "A blade sampling could not be loaded. Use method sampling."

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
        maxlen = 2*math.sqrt(sum(e**2))+0.03
        self._points = zeros((0,6))
        for side in sides:
            ex = math.sqrt(sum(side[6:9]**2))
            ey = math.sqrt(sum(side[9:12]**2))
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
        
    def make_model(self):
        try:
            makedirs('./Blade/'+self._model.model_type)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                print "BladeModeling::make_model - Problem creating Blade/mode_type folder"
                raise
        self._model._points = self._points
        w, RBFpoints, kernel = self._model.make()
        savez_compressed('Blade/'+self._model.model_type+'/'+self._name+'_'+kernel+'_w.npz', array=w)
        savez_compressed('Blade/'+self._model.model_type+'/'+self._name+'_'+kernel+'_points.npz', array=RBFpoints)
            

    def generate_trajectory(self, iter_surface):
        """ Method generate the coating trajectories. The trajectories are
        the intersection between two surfaces: the blade model, and the surface
        to be iterated, e.g. spheres (radius 1.425 to 3.770). The algorithm
        follows the marching method, documentation available in:
        http://www.mathematik.tu-darmstadt.de/~ehartmann/cdgen0104.pdf
        page 94 intersection surface - surface.
        """

        print "BladeModeling::generate_trajectory - Warning: this is a data-intensive computing and might freeze your computer."
        try:
            makedirs('./Blade/Trajectory')
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                print "BladeModeling::generate_trajectory - Problem creating Blade/Trajectory folder"
                raise
        try:
            Y = load('Blade/Trajectory/'+self._name+'.npz')
            Y = Y['array']
            tempY = []
            for y in Y:
                tempY.append(list(y))   
            Y = tempY
            Pd = Y[-1][-1]
            Rn = iter_surface.find_iter(Pd)
            Pd = mathtools.curvepoint(self._model, iter_surface, [Pd[0],Pd[1],Pd[2]])
            norm = self._model.df(Pd)
            norm /= sqrt(dot(norm,norm))
            Pd = [Pd[0],Pd[1],Pd[2],norm[0],norm[1],norm[2]]
            print "BladeModeling::generate_trajectory - Trajectories are loaded."
        except:
            print "BladeModeling::generate_trajectory - Trajectories could not be loaded."
            Y=[[]]
            Pd = initialPoint()
