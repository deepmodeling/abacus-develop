list(APPEND kernels_srcs
    blas_op.cpp
    lapack_op.cpp
    memory_op.cpp
)

if(ENABLE_CUDA_TOOLKIT)
    list(APPEND kernels_srcs
        cuda/memory_op.cu
        cuda/lapack_op.cu
        cuda/blas_op.cu
    )
endif()

if(ENABLE_HIP_TOOLKIT)
    hip_add_library(kernels_rocm STATIC
        rocm/memory_op.hip.cu
    )
endif()

add_library(kernels OBJECT ${kernels_srcs})