import numpy
a=numpy.load('grad_vepsl.npy')
b=numpy.load('OUT.autotest/deepks_stot.npy')
c=numpy.load('OUT.autotest/deepks_sbase.npy')
print(numpy.sum(numpy.absolute(a))+numpy.sum(numpy.absolute(b))+numpy.sum(numpy.absolute(c)))
