import coating
from openravepy import *
from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from math import *

# Get poinst for coating
nearlist = numpy.load('nearPointsByNumberOfPoints0_8.npz')
nearlist  = nearlist['array']

allangles = numpy.load('allangles0_8.npz')
allangles  = allangles['array']

allcouplesangles, alltriopoints = coating.pairingAllAngles(allangles,nearlist)
deltasT = coating.deltaTCalc(alltriopoints)

savez_compressed('alltriopoints0.npz', array=alltriopoints)
savez_compressed('deltasT0.npz', array=deltasT)
