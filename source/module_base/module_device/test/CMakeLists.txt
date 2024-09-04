AddTest(
  TARGET Module_Device_UTs
  LIBS ${math_libs} base device
  SOURCES memory_test.cpp device_test.cpp
)