import openravepy
from numpy import *
from openravepy import *
env=Environment()
env.SetViewer('qtcoin')
env.Load("/home/flash/Documents/openrave/kuka/arm_example.xml")
#env.Load("robots/kuka-kr30l16.zae")
robot = env.GetRobots()[0]
RaveSetDebugLevel(DebugLevel.Debug)

with env:
    ikmodel = databases.inversekinematics.InverseKinematicsModel(robot=robot,iktype=IkParameterization.Type.Transform6D)
    if not ikmodel.load():
        ikmodel.autogenerate()
    basemanip = interfaces.BaseManipulation(robot)
    reachmodel = databases.kinematicreachability.ReachabilityModel(robot)
    if not reachmodel.load():
        print "nao achou o reach"

    T_A = array([[ 1,  0,  0,   100], [  0,1,   0, 1.5], [  0,   0  , 1,   0], [  0,0,0,1]])
    initial_angles = [ 0,  0,  0,  0, 0,  0]
    robot.SetDOFValues(initial_angles,ikmodel.manip.GetArmIndices())


    iksol = ikmodel.manip.FindIKSolution(T_A,IkFilterOptions.CheckEnvCollisions)
    if iksol is None:
        print "position not reacheable"
    else:
        print "position reached" 
        robot.SetDOFValues(iksol,ikmodel.manip.GetArmIndices())

        

