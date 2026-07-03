abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)

if(USE_CUDA OR USE_ROCM)
    AddTest(
        TARGET MODULE_HSOLVER_KERNELS_Unittests
        LIBS parameter ${math_libs} base device
        SOURCES math_hegvd_test.cpp
    )
endif()

if(ENABLE_GOOGLEBENCH)
    AddTest(
        TARGET PERF_MODULE_HSOLVER_KERNELS
        LIBS parameter ${math_libs} base device
        SOURCES perf_math_kernel.cpp
    )
endif()