import numpy
a=numpy.load('OUT.autotest/deepks_otot.npy')
b=numpy.load('OUT.autotest/deepks_obase.npy')
print((a-b)[0][0])
