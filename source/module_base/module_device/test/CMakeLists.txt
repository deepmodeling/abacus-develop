AddTest(
  TARGET Module_Device_UTs
  LIBS parameter ${math_libs} base device 
  SOURCES memory_test.cpp device_test.cpp ../../tool_quit.cpp
)