list(APPEND container_srcs
    tensor.cpp
    tensor_map.cpp
    cpu_allocator.cpp
    tensor_buffer.cpp
    tensor_shape.cpp
    tensor_types.cpp
)

if(USE_CUDA)
    list(APPEND container_srcs gpu_allocator.cpp)
endif()

add_library(container OBJECT ${container_srcs})

add_subdirectory(kernels)

target_link_libraries(container kernels)

if(ENABLE_COVERAGE)
  add_coverage(container)
endif()

if(BUILD_TESTING)
  add_subdirectory(tests)
endif()
