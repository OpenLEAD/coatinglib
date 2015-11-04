import coating
import numpy
coatedarray = numpy.load('coatedpoints1_fullhd.npz')
coatedarray  = coatedarray['array']
nearlist = coating.nearPointsByNumberOfPoints(coatedarray)
numpy.savez_compressed('nearPointsByNumberOfPoints1_fullHD.npz', array=nearlist)
