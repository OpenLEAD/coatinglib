import coating
from openravepy import *
from numpy import *
from math import *

# Get poinst for coating
alltriopoints = load('alltriopoints0_HD.npz')
alltriopoints  = alltriopoints['array']

deltasT = load('deltasT0_HD.npz')
deltasT  = deltasT['array']

velocities,accelerations = coating.calculateLinearVelocitiesAndAccelerations(alltriopoints,deltasT)

savez_compressed('linearAcc0_HD.npz', array=accelerations)
savez_compressed('linearVel0_HD.npz', array=velocities)
