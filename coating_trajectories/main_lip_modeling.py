from rbf import RBF
import blade_modeling
from turbine import Turbine
from numpy import dot, sqrt, sum, concatenate, array, arange
import mathtools
from turbine_config import TurbineConfig, ConfigFileError
from visualizer import Visualizer
from os.path import join, exists, realpath
from os import makedirs, environ
from openravepy import matrixFromAxisAngle
from math import pi
from copy import deepcopy

def load_samples():
    blade.load_samples(join(name,'samples','samples.xml'))
    return

def save_samples():
    blade.save_samples(join(name,'samples'), name)
    return

def sampling_blade():
    blade.samples_delta = 0.005
    blade.sampling()

    # Extracting Lip
    Rs = sqrt(sum(blade.points[:,0:3]*blade.points[:,0:3],1))
    maxRs = 3.7560272702070847
    lip = blade.points[Rs>maxRs]
    lip_face = blade.points[(dot(blade.points[:,3:6],[1,0,0])>0.5)&(Rs>3.775)]
    lip_up = blade.points[(dot(blade.points[:,3:6],[0,0,1])>0.5)]
    blade.points = lip

    #Filtering
    blade.points = mathtools.filter_by_distance(blade.points, r = 0.009)
    
    return 

def plot_samples():
    load_samples()
    vis.plot(blade.points, 'sampling', ((1,0,0)))
    return

# ------------------------------------------------------------------------------------------------
# RBF Model
# ------------------------------------------------------------------------------------------------
def load_model():
    blade.load_model(join(name,'model','model.xml'))
    return

def save_model():
    blade.save_model(join(name,'model'), name)
    return 

def make_model():
    blade.make_model(rbf_model)
    return
        
# ------------------------------------------------------------------------------------------------
# Generate Trajectories
# ------------------------------------------------------------------------------------------------
def load_trajectories():
    blade.load_trajectory(join(name,'trajectory','trajectory.xml'))
    return


def save_trajectories():
    blade.save_trajectory(join(name,'model','model.xml'),
                          join(name,'trajectory'), name)
    return

def plot_trajectories():
    load_trajectories()
    for i in range(0,len(blade.trajectories),5):
        vis.plot(blade.trajectories[i], 'trajectories', ((0,0,1)))
    return

def generate_trajectories():

    load_samples()
    load_model()
    #Rs = sqrt(sum(blade.points[:,0:3]*blade.points[:,0:3],1))
    #lip_face = blade.points[(dot(blade.points[:,3:6],[1,0,0])>0.5)&(Rs>3.775)]
    lip_up = blade.points[(dot(blade.points[:,3:6],[0,0,1])>0.5)]

    #blade.points = lip_up
    #blade.make_model(rbf_model)

    class IterRBF(mathtools.IterSurface):
        def __init__(self, Rn0=0, stopR=0.05, coatingstep = 0.003, points = []):
            mathtools.IterSurface.__init__(self, Rn0, stopR, coatingstep)
            self.points = points
            self.rbf = RBF('r3', self.points)
            self.rbf.make()

        def update(self):
            self.rbf._points[:,0:3] = self.rbf._points[:,0:3] - self.coatingstep*self.rbf._points[:,3:6] 
            return

        def f(self,p):
            return self.rbf.f(p)

        def df(self,p):
            return self.rbf.df(p)

        def f_array(self,p):
            result = []
            for pi in p:
                result.append(self.f(pi))
            return result

        def find_iter(self, p0):
            n = abs(self.f(p0))
            if n<=1e-2:
                return
            
            for i in arange(self._Rn0, self.stopR, self.coatingstep):
                self.rbf._points[:,0:3] = self.rbf._points[:,0:3] - self.coatingstep*self.rbf._points[:,3:6]
                n = abs(self.f(p0))
                if n<=1e-2:
                    self._Rn = i+self.coatingstep
                    break
            else: raise ValueError('Outside rbf')
            return

        def criteria(self):
            return self._Rn<self.stopR

        def findnextparallel(self, p):
            self.find_iter(p)
            self.update()
            return

        def name(self):
            return 'IterRBF'

    #rbf_face = IterRBF(points=lip_face)
    p = blade.points[abs(blade.points[:,0])==min(abs(blade.points[:,0]))][0]
    normal_plane = array(p[3:6])
    plane = mathtools.Plane(max(lip_up[:,1]), min(lip_up[:,1]), k=1)#normal_plane=normal_plane)#k=0)
    #plane.findnextparallel(p)
    
    try:
        blade.generate_trajectories(plane)
    except KeyboardInterrupt:
        save_trajectories()
        raise

def reorganize_trajectories():
    load_trajectories()
    maxRs = 3.7560272702070847
    for i in range(0,len(blade.trajectories)):
        blade.trajectories[i] = array(blade.trajectories[i])
        Rs = sqrt(sum(blade.trajectories[i][:,0:3]*blade.trajectories[i][:,0:3],1))
        blade.trajectories[i] = (blade.trajectories[i])[(
            (dot(blade.trajectories[i][:,3:6],[0,0,1])>0.5)&(Rs>maxRs))|(
            (dot(blade.trajectories[i][:,3:6],[1,0,0])>0.5)&(Rs>3.775))]

    new_blade_trajectories = mathtools.trajectory_verticalization(
        blade.trajectories, step = 3)
    blade.trajectories = new_blade_trajectories
        
def filter_trajectories():
    load_trajectories()
    blade.trajectories = blade.filter_trajectory_opt()
    vis.plot_lists(blade.trajectories)
    return 

def save_filtered_trajectories():
    blade.save_trajectory(join(name,'model','model.xml'),
			  join(name+'_filtered','trajectory'),
			  name)
    return

def load_filtered_trajectories():
    blade.load_trajectory(join(name+'_filtered','trajectory','trajectory.xml'))
    return

def rotate_blade(T, directory):
    load_trajectories()
    blade.rotate_models(T)
    blade.trajectories = mathtools.rotate_trajectories(blade.trajectories, T)
    blade.save_model(join(directory,'model'), directory)
    blade.save_trajectory(join(directory,'model','model.xml'),
                          join(directory,'trajectory'), directory)
    return
    
if __name__ == "__main__":
    
    name = 'lip'#_filtered'
    dir_test = join(realpath('.'),'test')
    environ['OPENRAVE_DATA'] = str(dir_test)
    cfg = TurbineConfig.load('turbine_unittest.cfg','test')
    turb = Turbine(cfg)
    turb.env.Remove(turb.primary)
    turb.env.Remove(turb.secondary)

    rbf_model = RBF('r3')
    blade = blade_modeling.BladeModeling(turb, turb.blades[1])

    # Visualizer
    vis = Visualizer(turb.env)
    #plot_samples()
    #load_trajectories()
    #plot_trajectories()



    #sampling_blade()
    #save_samples()
    
    #make_model()
    #save_model()

    #generate_trajectories()
    #save_trajectories()

    filter_trajectories()
    #save_filtered_trajectories()


    
