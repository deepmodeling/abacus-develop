if(ENABLE_LCAO AND USE_CUDA)
  AddTest(
  TARGET gint_gpu_test
  LIBS parameter ${math_libs} psi base device
  SOURCES test_sph.cu test_sph.cpp
)
endif()