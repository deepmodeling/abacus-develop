remove_definitions(-D__MPI)

# Disable this test on Sunway (USE_SW) due to LAPACK symbol issues.
if(USE_SW)
  message(WARNING "Skip MODULE_BASE_CONTAINER_Unittests on USE_SW.")
else()
  AddTest(
    TARGET MODULE_BASE_CONTAINER_Unittests
    LIBS parameter container base device ${math_libs}
    SOURCES tensor_test.cpp tensor_shape_test.cpp allocator_test.cpp
    tensor_buffer_test.cpp tensor_map_test.cpp
    tensor_utils_test.cpp tensor_accessor_test.cpp
  )
endif()
