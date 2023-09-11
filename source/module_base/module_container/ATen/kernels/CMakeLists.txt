file(GLOB ATen_KERNELS_SRCS "*.cpp")

if(USE_CUDA)
    file(GLOB ATen_KERNELS_CUDA_SRCS "cuda/*.cu")
    set(ATen_CUDA_SRCS ${ATen_CUDA_SRCS} ${ATen_KERNELS_CUDA_SRCS} PARENT_SCOPE)
endif()

if(USE_ROCM)
    file(GLOB ATen_KERNELS_ROCM_SRCS "cuda/*.hip.cu")
    set(ATen_ROCM_SRCS ${ATen_ROCM_SRCS} ${ATen_KERNELS_ROCM_SRCS} PARENT_SCOPE)
endif()

set(ATen_CPU_SRCS ${ATen_CPU_SRCS} ${ATen_KERNELS_SRCS} PARENT_SCOPE)

if(BUILD_TESTING)
  add_subdirectory(test)
endif()