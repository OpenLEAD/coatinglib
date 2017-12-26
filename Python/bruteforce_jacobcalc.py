import coating
from openravepy import *
from numpy import *
from math import *
import time

env=Environment()
env.Load("/home/renan/workspace/coatinglib/motoman/motoman_mh122p.xml")
robot = env.GetRobots()[0]
manip = robot.GetActiveManipulator()

ikmodel = databases.inversekinematics.InverseKinematicsModel(robot=robot,iktype=IkParameterization.Type.Transform6D)
if not ikmodel.load():ikmodel.autogenerate()
        
armjoints = [j for j in robot.GetDependencyOrderedJoints() if j.GetJointIndex() in ikmodel.manip.GetArmIndices()]

Limits=[]
for joint in armjoints:
    limits = []
    limits.append(joint.GetLimits()[0][0])
    limits.append(joint.GetLimits()[1][0])
    if limits[0] == -limits[1]:
        limits[1]=0
    if limits[0]==-2*pi and limits[1]==0:
        limits[0]=-pi
    Limits.append(limits)

Limits = array(Limits)*180/pi
step=20
Limits=Limits/step

angles = []
for limit in Limits:
    lista=[]
    for i in range(int(limit[0]),int(limit[1])):
        lista.append(i*step*pi/180)
    angles.append(array(lista))

Ts=[]
Js=[]
Ss=[]
S_time = time.time() 
for S in angles[0]:
    S_time = time.time()
    for L in angles[1]:
        L_time = time.time()
        for U in angles[2]:
            U_time = time.time()
            for R in angles[3]:
                #R_time = time.time()
                for B in angles[4]:
                    #B_time = time.time()
                    for T in angles[5]:
                        robot.SetDOFValues([float(S),float(L),float(U),float(R),float(B),float(T)],ikmodel.manip.GetArmIndices())
                        Jpos = manip.CalculateJacobian()
                        Jori = manip.CalculateAngularVelocityJacobian()
                        J = concatenate((Jpos,Jori))
                        Ts.append(manip.GetTransform())
                        Js.append(J)
                        _, s, _ = linalg.svd(J, full_matrices=True)
                        Ss.append(s)
                    #print 'B_time = '+str(time.time()-B_time)
                #print 'R_time = '+str(time.time()-R_time)
            #tempo = time.time()-U_time            
            #print 'U_time = '+str(time.time()-U_time)
            #print 'Estimated time = '+str(tempo*len(angles[1])*len(angles[0]))
        print 'L_time = '+str(time.time()-L_time)
    print 'S_time = '+str(time.time()-S_time)
    
