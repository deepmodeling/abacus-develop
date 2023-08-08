import numpy
a=numpy.load('h_tot.npy')
b=numpy.load('h_base.npy')
print(numpy.sum(a-b))
