import numpy
a=numpy.load('OUT.autotest/deepks_gradvx.npy')
b=numpy.load('OUT.autotest/deepks_ftot.npy')
c=numpy.load('OUT.autotest/deepks_fbase.npy')
print(numpy.sum(numpy.absolute(a))+numpy.sum(numpy.absolute(b))+numpy.sum(numpy.absolute(c)))
