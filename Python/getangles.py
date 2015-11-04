import coating
import numpy

# Get poinst for coating
nearlist = numpy.load('nearPointsByNumberOfPoints0_fullHD.npz')
nearlist  = nearlist['array']

allangles = coating.computeAllAngularDistances(nearlist)
numpy.savez_compressed('allangles0_fullHD.npz', array=allangles)
