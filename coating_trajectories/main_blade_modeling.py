from rbf import RBF
import blade_modeling
from turbine import Turbine
from numpy import dot, sqrt, sum, concatenate, pi
import mathtools
from turbine_config import TurbineConfig
from visualizer import Visualizer
from os.path import join, realpath
from os import environ
from openravepy import matrixFromAxisAngle

def load_samples():
    blade.load_samples(join(name,'samples','samples.xml'))
    return

def save_samples():
    blade.save_samples(join(name,'samples'), name)
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
    #blade.points = mathtools.filter_by_distance(blade.points, 0.02)
    #lip = mathtools.filter_by_distance(lip, 0.01)
    #blade.points = concatenate((blade.points, lip))
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
    blade.save_trajectory(join(name,'model','model.xml'),
                          '', name)
    return

def plot_trajectories():
    load_trajectories()
    for i in range(0,len(blade.trajectories),5):
        vis.plot(blade.trajectories[i], 'trajectories', (0,0,1),1)
    return

def generate_trajectories():
    Rs = sqrt(sum(blade.points[:,0:3]*blade.points[:,0:3],1))
    sphere = mathtools.Sphere(max(Rs), min(Rs), 0.01)
    blade.generate_trajectories(sphere, vis)
##
##    try:
##        blade.generate_trajectories(sphere)
##    except KeyboardInterrupt:
##        save_trajectories()
##        raise

def filter_trajectories():
    load_trajectories()
    blade.trajectories = blade.filter_trajectory_opt()
    return 

def save_filtered_trajectories(trajectories):
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
    
    name = 'FACE/jiraublade_filtered'
    dir_test = join(realpath('.'),'test')
    environ['OPENRAVE_DATA'] = str(dir_test)
    cfg = TurbineConfig.load('turbine_unittest.cfg','test')
    turb = Turbine(cfg)
    
    turb.env.Remove(turb.primary)
    turb.env.Remove(turb.secondary)
    turb.env.Remove(turb.iris)
    turb.env.Remove(turb.rotor)
    turb.env.Remove(turb.runner_area)
    turb.env.Remove(turb.blades[0])
    turb.env.Remove(turb.blades[2])
    turb.env.Remove(turb.blades[1])
    turb.env.Remove(turb.robot)

    rbf_model = RBF('r3')
    blade = blade_modeling.BladeModeling(turb, turb.blades[3])

    # Visualizer
    vis = Visualizer(turb.env)
    #plot_samples()
    #plot_trajectories()

    sampling_blade()
    #load_samples()
    #save_samples()
    
    #make_model()
    #load_model()
    #save_model()
    #T = matrixFromAxisAngle([pi, 0, 0])
    #blade.rotate_models(T)

    #x = raw_input('wait')
    #generate_trajectories()
    #save_trajectories()

    #filter_trajectories()
    #save_filtered_trajectories()

    #rotate_blade(matrixFromAxisAngle([-pi/4, 0, 0]), 'jiraublade_hd_-45')

    
    
