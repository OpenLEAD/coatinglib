#!/usr/bin/env python
import time
import openravepy
from numpy import *
from openravepy import *

env = Environment() # create openrave environment
env.SetViewer('qtcoin') # attach viewer (optional)
env.Load('kuka_test.env.xml') # load a simple scene
robot = env.GetRobots()[0]
robot2 = env.GetRobots()[1]
robot3 = env.GetRobots()[2]
time.sleep(0.1) # give time for environment to update

with env:
    ikmodel = databases.inversekinematics.InverseKinematicsModel(robot=robot,iktype=IkParameterization.Type.Transform6D)
    if not ikmodel.load():
        ikmodel.autogenerate()


    reachmodel = databases.kinematicreachability.ReachabilityModel(robot)
    print reachmodel.getfilename()
    if not reachmodel.load():
        print "nao achou o reach"

    #reachmodel = databases.kinematicreachability.ReachabilityModel(robot)

    basemanip = interfaces.BaseManipulation(robot)
    taskmanip = interfaces.TaskManipulation(robot)

    #initial_angles = [ 0,  0,  -0.55,  0, 0,  0]
    #obot.SetDOFValues(initial_angles,ikmodel.manip.GetArmIndices())
    #robot.WaitForController(0)

    #each limit of the blade, begining from the top left corner and going clockwise (with 20cm distance)
    T_A = array([[ -1,  0,  0,   1.35], [  0,1,   0, 0.7], [  0,   0  , 1,   -1.31], [  0,0,0,1]])
    T_B = array([[ 1,  0,  0,   1.30], [  0,1,   0, 1.6], [  0,   0  , 1,   -3.25], [  0,0,0,1]])
    T_C = array([[ 1,  0,  0,   2.2], [  0,1,   0, 0], [  0,   0  , 1,   2.4], [  0,0,0,1]])
    T_D = array([[ 1,  0,  0,   2.5], [  0,1,   0, -0.7], [  0,   0  , 1,   2.5], [  0,0,0,1]])
    T_E = array([[ 1,  0,  0,   2.5], [  0,1,   0, -1.1], [  0,   0  , 1,   2.8], [  0,0,0,1]])
    T_F = array([[ 1,  0,  0,   1.6], [  0,1,   0, 1.7], [  0,   0  , 1,   0.6], [  0,0,0,1]])
    T_G = array([[ 1,  0,  0,   2.1], [  0,1,   0, 0], [  0,   0  , 1,   0.1], [  0,0,0,1]])
    T_H = array([[ 1,  0,  0,   2.5], [  0,1,   0, -1.5], [  0,   0  , 1,   0.4], [  0,0,0,1]])

    corners = [T_A,T_B,T_C,T_D,T_E,T_F,T_G,T_H]
    
    #time.sleep(10)

    #basemanip.MoveToHandPosition(matrices=[T_A],seedik=10)
    #robot.WaitForController(0)




    #for corner in corners:
       # basemanip.MoveToHandPosition(matrices=[corner],seedik=10)
     #   #robot.WaitForController(0)


    iksol = ikmodel.manip.FindIKSolution(T_A,IkFilterOptions.CheckEnvCollisions)
    if iksol is None:
        print "position not reacheable"
    else:
        print "position reached" 
        robot.SetDOFValues(iksol,ikmodel.manip.GetArmIndices())


    #for i in corners:
    #    iksol = ikmodel.manip.FindIKSolution(i,IkFilterOptions.CheckEnvCollisions)
    #    if iksol is None:
    #        print "position not reacheable"
    #    else:
    #        print "position reached" 
    #        basemanip.MoveToHandPosition(matrices=[i],seedik=10)
    #        #robot.SetDOFValues(iksol,ikmodel.manip.GetArmIndices())
    #        robot.WaitForController(0)
    
    #res = basemanip.MoveToHandPosition(matrices=[Tstart],seedik=10)
    

    #time.sleep(3) # give time for environment to update
    #Tgoal = array([[ 1,  0,  0,   2], [  0,1,   0, -2], [  0,   0  , 1,   2], [  0,0,0,1]])
    #sol2 = ikmodel.manip.FindIKSolution(Tgoal,IkFilterOptions.CheckEnvCollisions)
    #robot.SetDOFValues(sol2,ikmodel.manip.GetArmIndices())

#self = openravepy.examples.constraintplanning.ConstraintPlanning(robot)
