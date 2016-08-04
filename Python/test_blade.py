import RBF_c
import blade_c
import math
import mathtools
from turbine import Turbine


rbf = RBF_c.RBF('jiraublade','r3')

turbine = Turbine('turbine_std.cfg', False)

blademodel = blade_c.BladeModeling('jiraublade', rbf, turbine, True)

blademodel.sampling()
##
blademodel.make_model()

sphere = mathtools.sphere(turbine.model.runner_radius,
                          turbine.model.nose_radius,
                          turbine.coating.parallel_gap
                          )
blademodel.generate_trajectory(sphere)
