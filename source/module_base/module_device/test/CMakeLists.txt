AddTest(
  TARGET BASE_MODULE_DEVICE_UnitTests
  LIBS parameter ${math_libs} base device 
  SOURCES memory_test.cpp device_test.cpp
)