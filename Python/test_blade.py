from openravepy import Environment
import RBF_c
import blade_c
import math

env=Environment()
env.Load("../Turbina/env_mh12_0_16.xml")


rbf = RBF_c.RBF('jiraublade','r3')

blade = blade_c.Blade('jiraublade', env)
blademodel = blade_c.BladeModeling('jiraublade', rbf, env, blade.blade)
blademodel.sampling(bladerotation=[-27*math.pi/180,'y'])
