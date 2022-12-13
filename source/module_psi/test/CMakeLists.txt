AddTest(
  TARGET Module_Psi_UTs
  LIBS ${math_libs} base device
  SOURCES memory_test.cpp device_test.cpp
          ../../src_parallel/parallel_reduce.cpp 
          ../../src_parallel/parallel_global.cpp 
)