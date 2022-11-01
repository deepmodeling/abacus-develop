if (USE_CUDA)
AddTest(
  TARGET Module_Psi_UTs
  LIBS ${math_libs}
  SOURCES memory_test.cpp device_test.cpp
          ../src/memory_psi.cpp
          ../src/device.cpp
          ../src/cuda/memory.cu
)
else()
AddTest(
  TARGET Module_Psi_UTs
  LIBS ${math_libs}
  SOURCES memory_test.cpp device_test.cpp
          ../src/memory_psi.cpp
          ../src/device.cpp
)
endif()