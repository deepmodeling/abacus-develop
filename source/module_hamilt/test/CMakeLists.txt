AddTest(
  TARGET Hamilt_UTs
  LIBS ${math_libs} device base
  SOURCES ekinetic_test.cpp nonlocal_test.cpp veff_test.cpp meta_test.cpp
          ../../src_parallel/parallel_reduce.cpp
          ../../src_parallel/parallel_global.cpp
          ../src/ekinetic.cpp
          ../src/nonlocal.cpp
          ../src/veff.cpp
          ../src/cuda/ekinetic.cu
          ../src/cuda/nonlocal.cu
          ../src/cuda/veff.cu
)
