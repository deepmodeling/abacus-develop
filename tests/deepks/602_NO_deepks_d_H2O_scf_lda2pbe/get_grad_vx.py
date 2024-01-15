import numpy
a=numpy.load('grad_vx.npy')
b=numpy.load('f_tot.npy')
c=numpy.load('f_base.npy')
print(numpy.sum(numpy.absolute(a))+numpy.sum(numpy.absolute(b))+numpy.sum(numpy.absolute(c)))
