import numpy
a=numpy.load('OUT.autotest/deepks_htot.npy')
b=numpy.load('OUT.autotest/deepks_hbase.npy')
print(numpy.sum(a-b))
