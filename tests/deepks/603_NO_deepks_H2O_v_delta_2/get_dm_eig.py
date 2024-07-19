import numpy
a=numpy.load('dm_eig.npy')
b=numpy.load('e_tot.npy')
c=numpy.load('e_base.npy')
print(numpy.sum(numpy.absolute(a))+numpy.sum(b)+numpy.sum(c))
