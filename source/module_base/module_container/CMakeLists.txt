# FetchContent_MakeAvailable was not added until CMake 3.14
cmake_minimum_required(VERSION 3.14)

set(ATen_CPU_SRCS)
set(ATen_CUDA_SRCS)
set(ATen_ROCM_SRCS)
set(ATen_CPU_DEPENDENCY_LIBS)
set(ATen_CUDA_DEPENDENCY_LIBS)
set(ATen_ROCM_DEPENDENCY_LIBS)

add_subdirectory(base/core)
add_subdirectory(ATen/core)
add_subdirectory(ATen/kernels)

if(USE_ROCM)
  hip_add_library(container_rocm STATIC ${ATen_ROCM_SRCS})
  set(ATen_ROCM_DEPENDENCY_LIBS container_rocm)
endif()

add_library(container OBJECT ${ATen_CPU_SRCS} ${ATen_CUDA_SRCS})

target_link_libraries(container PUBLIC 
  ${ATen_CPU_DEPENDENCY_LIBS} ${ATen_CUDA_DEPENDENCY_LIBS} ${ATen_ROCM_DEPENDENCY_LIBS})

if(ENABLE_COVERAGE)
  add_coverage(container)
endif()

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
