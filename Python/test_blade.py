import RBF
import blade
import math
import mathtools
from turbine import Turbine

name = "jiraublade_t"

rbf = RBF.RBF(name,'r3')

turbine = Turbine('turbine_std.cfg', False)

blademodel = blade.BladeModeling(name, rbf, turbine, True)

blademodel.sampling()
##
#blademodel.make_model()

#sphere = mathtools.sphere(turbine.model.runner_radius,
#                          turbine.model.nose_radius,
#                          turbine.coating.parallel_gap
#                          )
#blademodel.generate_trajectory(sphere)
