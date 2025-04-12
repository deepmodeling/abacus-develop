import numpy
import sys
file=sys.argv[1]
a=numpy.load(file)
print(numpy.sum(a))
