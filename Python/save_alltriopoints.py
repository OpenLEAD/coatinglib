import coating
from openravepy import *
from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from math import *

# Get poinst for coating
nearlist = numpy.load('nearPointsByNumberOfPoints0_HD.npz')
nearlist  = nearlist['array']

allangles = numpy.load('allangles0_HD.npz')
allangles  = allangles['array']

allcouplesangles, alltriopoints = coating.pairingAllAngles(allangles,nearlist)
deltasT = coating.deltaTCalc(alltriopoints)

savez_compressed('alltriopoints0_HD.npz', array=alltriopoints)
savez_compressed('deltasT0_HD.npz', array=deltasT)
