remove_definitions(-D__DEEPKS)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)

install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

AddTest(
  TARGET io_input_test
  LIBS ${math_libs} base device
  SOURCES input_test.cpp ../input.cpp
)

AddTest(
  TARGET io_winput_test
  LIBS ${math_libs} base device
  SOURCES winput_test.cpp ../winput.cpp
)

AddTest(
  TARGET io_output_test
  LIBS ${math_libs} base device
  SOURCES output_test.cpp ../output.cpp
)

AddTest(
  TARGET io_rwstream_test
  SOURCES rwstream_test.cpp ../rwstream.cpp
)
