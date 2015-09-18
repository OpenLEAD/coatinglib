import coating
import numpy

# Get poinst for coating
nearlist = numpy.load('nearPointsByNumberOfPoints0_HD.npz')
nearlist  = nearlist['array']

allangles = coating.computeAllAngularDistances(nearlist)
numpy.savez_compressed('allangles0_HD.npz', array=allangles)
