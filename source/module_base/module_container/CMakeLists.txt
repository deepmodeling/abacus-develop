set(ATen_CPU_SRCS)
set(ATen_CUDA_SRCS)
set(ATen_ROCM_SRCS)
set(ATen_CPU_DEPENDENCY_LIBS)
set(ATen_CUDA_DEPENDENCY_LIBS)
set(ATen_ROCM_DEPENDENCY_LIBS)

add_subdirectory(base/core)
add_subdirectory(ATen/core)
add_subdirectory(ATen/kernels)

# Pass source, includes, and libs to parent
set(ATen_CPU_SRCS ${ATen_CPU_SRCS} PARENT_SCOPE)
set(ATen_CUDA_SRCS ${ATen_CUDA_SRCS} PARENT_SCOPE)
set(ATen_ROCM_SRCS ${ATen_ROCM_SRCS} PARENT_SCOPE)
set(ATen_CPU_DEPENDENCY_LIBS ${ATen_CPU_DEPENDENCY_LIBS} PARENT_SCOPE)
set(ATen_CUDA_DEPENDENCY_LIBS ${ATen_CUDA_DEPENDENCY_LIBS} PARENT_SCOPE)
set(ATen_ROCM_DEPENDENCY_LIBS ${ATen_ROCM_DEPENDENCY_LIBS} PARENT_SCOPE)

add_library(container OBJECT ${ATen_CPU_SRCS} ${ATen_CUDA_SRCS})

if(ENABLE_COVERAGE)
  add_coverage(container)
endif()

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
