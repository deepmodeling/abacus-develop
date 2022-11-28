if (USE_CUDA)
AddTest(
  TARGET Module_Psi_UTs
  LIBS ${math_libs} base
  SOURCES memory_test.cpp device_test.cpp
          ../src/memory_psi.cpp
          ../src/device.cpp
          ../src/cuda/memory.cu
          ../../src_parallel/parallel_reduce.cpp 
          ../../src_parallel/parallel_global.cpp 
)
else()
AddTest(
  TARGET Module_Psi_UTs
  LIBS ${math_libs} base
  SOURCES memory_test.cpp device_test.cpp
          ../src/memory_psi.cpp
          ../src/device.cpp
          ../../src_parallel/parallel_reduce.cpp 
          ../../src_parallel/parallel_global.cpp 
)
endif()