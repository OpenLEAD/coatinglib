from rbf import RBF
import blade_modeling
from turbine import Turbine
from numpy import dot, sqrt, sum, concatenate, array
import mathtools
from turbine_config import TurbineConfig, ConfigFileError
from visualizer import Visualizer
from os.path import join, exists
from os import makedirs, environ

def load_samples():
    blade.load_samples(name+'/samples/samples.xml')
    return

def save_samples():
    blade.save_samples(join(name,'samples'), name)
    return

def plot_samples():
    vis.plot(blade.points, 'sampling', ((1,0,0)))
    return

def sampling_blade():
    blade.sampling()

    # Removing points which normal vectors' directions are close to [0,1,0] and samples above y = -2
    y = [0,1,0]
    #blade.points = blade.points[~((dot(blade.points[:,3:6],y)>=0.7)&(blade.points[:,1]>-2))]

    # Removing points which normal vectors' directions are close to [0,-1,0] 
    #blade.points = blade.points[~(dot(blade.points[:,3:6],y)<=-0.7)]

    # Removing Lip
    Rs = sqrt(sum(blade.points[:,0:3]*blade.points[:,0:3],1))
    maxRs = max(Rs)

    lip = blade.points[(Rs>max(Rs)-0.1)&(dot(blade.points[:,3:6],y)>=0.3)]
    blade.points = blade.points[~((Rs>max(Rs)-0.1)&(dot(blade.points[:,3:6],y)>=0.3))]

    # --------------------------
    # Filtering
    # --------------------------
    # Moving through normal
    #blade.points[:,0:3] = blade.points[:,0:3]+0.2*blade.points[:,3:6]
    #blade.points[:,0:3] = 2*blade.points[:,0:3]

    blade.points = blade.filter_by_distance(blade.points, 0.02)
    lip = blade.filter_by_distance(lip, 0.01)

    #blade.points[:,0:3] = 0.5*blade.points[:,0:3]
    #blade.points[:,0:3] = blade.points[:,0:3]-0.2*blade.points[:,3:6]

    blade.points = concatenate((blade.points, lip))
    return 

# ------------------------------------------------------------------------------------------------
# RBF Model
# ------------------------------------------------------------------------------------------------
def load_model():
    blade.load_model(name+'/model/model.xml')
    return

def save_model():
    blade.save_model(join(name,'model'), name)
    return 

def make_model():
    s = mathtools.Sphere(0,0,1)
    blade.make_model(rbf_model, s)
    return
        
# ------------------------------------------------------------------------------------------------
# Generate Trajectories
# ------------------------------------------------------------------------------------------------
def load_trajectories():
    blade.load_trajectory(join(name,'trajectory','trajectory.xml'))
    return


def save_trajectories():
    blade.save_trajectory(join(name,'model/model.xml'),
                          join(name,'trajectory'), name)
    return

def plot_trajectories():
    for trajectory in blade.trajectories:
        vis.plot(trajectory, 'trajectories', ((0,0,1)))
    return

def generate_trajectories():
    Rs = sqrt(sum(blade.points[:,0:3]*blade.points[:,0:3],1))
    sphere = mathtools.Sphere(max(Rs), min(Rs), 0.003)

    try:
        blade.generate_trajectories(sphere)
    except KeyboardInterrupt:
        save_trajectories()
        raise

def filter_trajectories():
    return blade.filter_trajectory_opt()

def save_filtered_trajectories():
    blade.trajectories = trajectories
    blade.save_trajectory(join(name,'model','model.xml'),
			  join(name+'_filtered','trajectory'),
			  name)
    return

def load_filtered_trajectories():
    blade.load_trajectory(join(name+'_filtered','trajectory','trajectory.xml'))
    return

if __name__ == "__main__":
    
    name = 'jiraublade_hd_filtered'
    turbconf = TurbineConfig.load("turbine_unittest.cfg", 'test/')
    environ['OPENRAVE_DATA'] = '/home/renan/git/planning-coating_trajectories/coating_trajectories/test'
    turb = Turbine(turbconf)
    turb.env.Remove(turb.primary)
    turb.env.Remove(turb.secondary)

    rbf_model = RBF('r3')
    blade = blade_modeling.BladeModeling(turb, turb.blades[0])

    #sampling_blade()
    #save_samples()
    #load_samples()
    
    #make_model()
    #save_model()
    #load_model()

    #generate_trajectories()
    #save_trajectories()
    load_trajectories()

    #trajectories = filter_trajectories()
    #save_filtered_trajectories()
    #load_filtered_trajectories()

    # Visualizer
    vis = Visualizer(turb.env)
    #plot_samples()
    plot_trajectories()
    
