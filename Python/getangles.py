import coating
import numpy

# Get poinst for coating
nearlist = numpy.load('nearPointsByNumberOfPoints1_8.npz')
nearlist  = nearlist['array']

allangles = coating.computeAllAngularDistances(nearlist)
numpy.savez_compressed('allangles1_8.npz', array=allangles)
