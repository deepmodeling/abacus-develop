remove_definitions(-D__MPI)
remove_definitions(-D__LCAO)

install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

AddTest(
  TARGET esolver_dp_test
  LIBS ${math_libs} base_serial device
  SOURCES esolver_dp_test.cpp ../esolver_dp.cpp
)
