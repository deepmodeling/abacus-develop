remove_definitions(-D__MPI)

AddTest(
  TARGET BASE_MODULE_CONTAINER_container_uts
  LIBS parameter container base device ${math_libs}
  SOURCES tensor_test.cpp tensor_shape_test.cpp allocator_test.cpp
  tensor_buffer_test.cpp tensor_map_test.cpp
  tensor_utils_test.cpp tensor_accessor_test.cpp
)
