abacus_disable_feature_definitions(__MPI)

addtest(
  TARGET MODULE_BASE_KERNELS_Unittests
  LIBS parameter base device
  SOURCES test_math_ylm_op.cpp test_math_kernel.cpp)
