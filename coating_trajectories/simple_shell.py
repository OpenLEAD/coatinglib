from turbine import Turbine
from turbine_config import TurbineConfig, ConfigFileError
import os
from os.path import join, isfile, realpath

dir_test = join(realpath('.'),'test')
os.environ['OPENRAVE_DATA'] = str(dir_test)
cfg = TurbineConfig.load('turbine_unittest.cfg','test')
turb = Turbine(cfg)
##import db
##bases = db.get_bases('db_bases_to_num.pkl')

import visualizer
vis = visualizer.Visualizer(turb.env)
