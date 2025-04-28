import numpy
import sys
file_tot=sys.argv[1]
file_base=sys.argv[2]
a=numpy.load(file_tot)
b=numpy.load(file_base)
print(numpy.sum(numpy.absolute(a-b)))
