import os
import turbine_config
import turbine
import visualizer
import numpy as np
import rail_place
import planning
import blade_modeling


dir_test = os.path.join(os.path.realpath('.'),'test')
os.environ['OPENRAVE_DATA'] = str(dir_test)

cfg = turbine_config.TurbineConfig.load('turbine_unittest.cfg','test')
turb = turbine.Turbine(cfg)

vz = visualizer.Visualizer(turb.env)

RP = rail_place.rand_rail(turb,20)
rpt = rail_place.RailPlace((.80,0.15,0))

folder = "jiraublade"
blade = blade_modeling.BladeModeling(turb, turb.blades[0])

#blade.load_trajectory(os.path.join(folder,"trajectory/trajectory.xml"))

##for trajectory in blade.trajectories:
##    vz.plot(trajectory)

#position trajectories
Ro = np.transpose(np.array([[-1,0,0],[0,0,1],[0,1,0]]))
To = np.eye(4)
To[0:3,0:3]=Ro
x=None

DB_dict = dict()

for rp in RP:
    
    turb.place_rail(rp)
    turb.place_robot(rp)
    turb.robot.SetTransform(np.dot(turb.robot.GetTransform(),To))
    
    if turb.check_rail_collision():
        continue
    
    if turb.check_robot_collision():
        continue

    filtered_trajectories = path_filters.filter_trajectories(turb, blade.trajectories)

    counter = 0
    for filtered_trajectory in filtered_trajectories:
        for filtered_trajectory_part in filtered_trajectory:
            try:
                lower, _, _ = compute_first_feasible_point(turb, filtered_trajectory_part)
            except:
                continue
            
            upper = lower + len(planning.compute_robot_joints(turb, filtered_trajectory_part, lower))
            
            for point in filtered_trajectory_part[lower:upper]:
                DB_dict[tuple(point[0:3])] = (DB_dict.get(tuple(point[0:3]),set()) |
                                              set([tuple(rp.getPSAlpha())]))

        print counter
        counter += 1


        
