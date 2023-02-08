import numpy
a=numpy.load('o_tot.npy')
b=numpy.load('o_base.npy')
print((a-b)[0][0])
