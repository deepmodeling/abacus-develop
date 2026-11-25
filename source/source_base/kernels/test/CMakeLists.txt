abacus_disable_feature_definitions(__MPI)

addtest(
  TARGET MODULE_BASE_KERNELS_Unittests
  LIBS parameter base device
  SOURCES math_ylm_op_test.cpp math_kernel_test.cpp)
