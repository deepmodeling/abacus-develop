AddTest(
  TARGET Kernels_UTs
  LIBS ${math_libs} base device
  SOURCES memory_op_test.cpp device_test.cpp
)