if (USE_CUDA)
AddTest(
  TARGET Hamilt_UTs
  LIBS ${math_libs}
  SOURCES ekinetic_test.cpp nonlocal_test.cpp
          ../src/ekinetic.cpp
          ../src/nonlocal.cpp
          ../src/cuda/ekinetic.cu
          ../src/cuda/nonlocal.cu
          ../../module_psi/src/memory.cpp
          ../../module_psi/src/cuda/memory.cu
)
else()
AddTest(
  TARGET Hamilt_UTs
  LIBS ${math_libs}
  SOURCES ekinetic_test.cpp nonlocal_test.cpp
          ../src/ekinetic.cpp
          ../src/nonlocal.cpp
          ../../module_psi/src/memory.cpp
)
endif()