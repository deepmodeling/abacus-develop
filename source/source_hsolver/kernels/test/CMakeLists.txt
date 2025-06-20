remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)

if(USE_CUDA OR USE_ROCM)
    AddTest(
        TARGET MODULE_HSOLVER_KERNELS_Unittests
        LIBS parameter ${math_libs} base device
        SOURCES math_dngvd_test.cpp
    )
endif()

if(ENABLE_GOOGLEBENCH)
    AddTest(
        TARGET PERF_MODULE_HSOLVER_KERNELS
        LIBS parameter ${math_libs} base device
        SOURCES perf_math_kernel.cpp
    )
endif()