remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)

if(USE_CUDA OR USE_ROCM)
    AddTest(
        TARGET Hsolver_Kernels_UTs
        LIBS ${math_libs} base device
        SOURCES math_kernel_test.cpp math_dngvd_test.cpp
        ../../../src_parallel/parallel_global.cpp ../../../src_parallel/parallel_common.cpp
        ../../../src_parallel/parallel_reduce.cpp
    )
elseif()
    AddTest(
        TARGET Hsolver_Kernels_UTs
        LIBS ${math_libs} base device
        SOURCES math_kernel_test.cpp ../../../src_parallel/parallel_global.cpp
        ../../../src_parallel/parallel_common.cpp  ../../../src_parallel/parallel_reduce.cpp
    )
endif()