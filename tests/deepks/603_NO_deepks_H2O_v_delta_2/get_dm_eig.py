import numpy
a=numpy.load('OUT.autotest/deepks_dm_eig.npy')
b=numpy.load('OUT.autotest/deepks_etot.npy')
c=numpy.load('OUT.autotest/deepks_ebase.npy')
print(numpy.sum(numpy.absolute(a))+numpy.sum(b)+numpy.sum(c))
