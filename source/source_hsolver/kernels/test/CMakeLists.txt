abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)

if(USE_CUDA OR USE_ROCM)
  addtest(
    TARGET MODULE_HSOLVER_KERNELS_Unittests
    LIBS parameter base device
    SOURCES math_hegvd_test.cpp)
endif()

if(ENABLE_GOOGLEBENCH)
  addtest(
    TARGET PERF_MODULE_HSOLVER_KERNELS
    LIBS parameter base device
    SOURCES perf_math_kernel.cpp)
endif()
