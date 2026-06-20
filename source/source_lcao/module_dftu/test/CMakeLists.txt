remove_definitions(-D__CUDA)

AddTest(
  TARGET dftu_pw_test
  LIBS ${math_libs} base device parameter
  SOURCES dftu_pw_test.cpp
)

AddTest(
  TARGET dftu_core_test
  LIBS ${math_libs} base device
  SOURCES dftu_core_test.cpp
)

AddTest(
  TARGET dftu_operator_test
  LIBS ${math_libs} base device
  SOURCES dftu_operator_test.cpp
)
