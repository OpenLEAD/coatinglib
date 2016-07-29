import RBF_c
import blade_c
import math
import mathtools
from turbine import Turbine


rbf = RBF_c.RBF('jiraublade','r3')

turbine = Turbine('turbine_std.cfg', False)

blademodel = blade_c.BladeModeling('jiraublade', rbf, turbine, False)

##blademodel.sampling()
##
#blademodel.make_model()

sphere = mathtools.sphere(turbine.environment.rotor_radius, turbine.environment.nose_radius,
                          turbine.coating.parallel_gap)
blademodel.generate_trajectory(sphere)
