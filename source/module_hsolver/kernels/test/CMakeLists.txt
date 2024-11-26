remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)

if(USE_CUDA OR USE_ROCM)
    AddTest(
        TARGET Hsolver_Kernels_UTs
        LIBS parameter ${math_libs} base device
        SOURCES math_kernel_test.cpp math_dngvd_test.cpp
    )
elseif()
    AddTest(
        TARGET Hsolver_Kernels_UTs
        LIBS parameter ${math_libs} base device
        SOURCES math_kernel_test.cpp ../../../module_base/blas_connector.cpp
    )
endif()

if(ENABLE_GOOGLEBENCH)
    AddTest(
        TARGET Perf_Hsolver_Kernels
        LIBS parameter ${math_libs} base device
        SOURCES perf_math_kernel.cpp
    )
endif()