import coating
import numpy
coatedarray = numpy.load('coatedpoints0.npz')
coatedarray  = coatedarray['array']
nearlist = coating.nearPointsByNumberOfPoints(coatedarray)
numpy.savez_compressed('nearPointsByNumberOfPoints0_HD.npz', array=nearlist)
