import numpy
a=numpy.load('grad_vepsl.npy')
b=numpy.load('s_tot.npy')
c=numpy.load('s_base.npy')
print(numpy.sum(numpy.absolute(a))+numpy.sum(numpy.absolute(b))+numpy.sum(numpy.absolute(c)))
