if (USE_CUDA)
AddTest(
  TARGET Module_Psi_UTs
  LIBS ${math_libs}
  SOURCES memory_test.cpp device_test.cpp
          ../../module_psi/src/memory_psi.cpp
          ../../module_psi/src/device.cpp
          ../../module_psi/src/cuda/memory.cu
)
else()
AddTest(
  TARGET Module_Psi_UTs
  LIBS ${math_libs}
  SOURCES memory_test.cpp device_test.cpp
          ../../module_psi/src/memory_psi.cpp
          ../../module_psi/src/device.cpp
)
endif()