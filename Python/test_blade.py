from openravepy import Environment
import RBF_c
import blade_c
import math
import mathtools


env=Environment()
env.Load("../Turbina/env_mh12_0_16.xml")
env.SetViewer('qtcoin')

rbf = RBF_c.RBF('jiraublade','r3')

blade = blade_c.Blade('jiraublade', env, 'pa1')
blademodel = blade_c.BladeModeling('jiraublade', rbf, env, blade.blade)

##blademodel.sampling(bladerotation=[-27*math.pi/180,'y'])
##
#blademodel.make_model()

sphere = mathtools.sphere()
blademodel.generate_trajectory(sphere)
