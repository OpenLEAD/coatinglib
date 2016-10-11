import os
dir_test = os.path.join(os.path.realpath('.'),'test')
os.environ['OPENRAVE_DATA'] = str(dir_test)
import turbine_config
cfg = turbine_config.TurbineConfig.load('turbine_unittest.cfg','test')
import turbine
turb = turbine.Turbine(cfg)
import visualizer
vz = visualizer.Visualizer(turb.env)
import numpy as np
import rail_place
RP = rail_place.rand_rail(turb,20)
rpt = rail_place.RailPlace((.80,0.15,0))

#position trajectories

for rp in RP:

    turb.place_rail(rp)
    if turb.check_rail_collision():
        continue
    
    turb.place_robot(rp)
    if turb.check_robot_collision():
        continue

    filtered_trajectories = filter_trajectories(turb,trajectories)
    
# filter PCL
# follow parallels (verifing limits) and saving points
# save railplace info to points

#loop

