AddTest(
  TARGET MODULE_BASE_DEVICE_Unittests
  LIBS parameter ${math_libs} base device 
  SOURCES memory_test.cpp device_test.cpp
)