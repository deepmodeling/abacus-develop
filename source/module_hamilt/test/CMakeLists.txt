if (USE_CUDA)
AddTest(
  TARGET Hamilt_UTs
  LIBS ${math_libs}
  SOURCES ekinetic_test.cpp nonlocal_test.cpp veff_test.cpp
          ../src/ekinetic.cpp
          ../src/nonlocal.cpp
          ../src/veff.cpp
          ../src/cuda/ekinetic.cu
          ../src/cuda/nonlocal.cu
          ../src/cuda/veff.cu
          ../../module_psi/src/memory_psi.cpp
          ../../module_psi/src/cuda/memory.cu
)
else()
AddTest(
  TARGET Hamilt_UTs
  LIBS ${math_libs}
  SOURCES ekinetic_test.cpp nonlocal_test.cpp veff_test.cpp
          ../src/ekinetic.cpp
          ../src/nonlocal.cpp
          ../src/veff.cpp
          ../../module_psi/src/memory_psi.cpp
)
endif()