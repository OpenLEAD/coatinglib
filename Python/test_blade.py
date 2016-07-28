import RBF_c
import blade_c
import math
import mathtools
from turbine import Turbine


rbf = RBF_c.RBF('jiraublade','r3')

turbine = Turbine('turbine_std.cfg', False)

blademodel = blade_c.BladeModeling('jiraublade', rbf, turbine, False)

##blademodel.sampling(bladerotation=[-27*math.pi/180,'y'])
##
#blademodel.make_model()

sphere = mathtools.sphere()
blademodel.generate_trajectory(sphere, True)
