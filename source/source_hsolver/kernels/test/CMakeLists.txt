remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)

if(USE_CUDA OR USE_ROCM)
  addtest(
    TARGET
    MODULE_HSOLVER_KERNELS_Unittests
    LIBS
    parameter
    ${math_libs}
    base
    device
    SOURCES
    math_hegvd_test.cpp)
endif()

if(ENABLE_GOOGLEBENCH)
  addtest(
    TARGET
    PERF_MODULE_HSOLVER_KERNELS
    LIBS
    parameter
    ${math_libs}
    base
    device
    SOURCES
    perf_math_kernel.cpp)
endif()
