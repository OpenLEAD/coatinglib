import coating
import numpy
coatedarray = numpy.load('coatedpoints.npz')
coatedarray  = coatedarray['array']
nearlist = coating.nearPointsByNumberOfPoints(coatedarray[1])
numpy.savez_compressed('nearPointsByNumberOfPoints1_8.npz', array=nearlist)
