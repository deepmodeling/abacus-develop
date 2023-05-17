remove_definitions(-D__MPI)
remove_definitions(-D__LCAO)

AddTest(
  TARGET esolver_dp_test
  LIBS ${math_libs} base device
  SOURCES esolver_dp_test.cpp ../esolver_dp.cpp
)