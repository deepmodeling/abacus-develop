remove_definitions(-D__MPI)
remove_definitions(-D__LCAO)

install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

AddTest(
  TARGET esolver_dp_test
  LIBS parameter ${math_libs} base device 
  SOURCES esolver_dp_test.cpp ../esolver_dp.cpp ../../module_io/cif_io.cpp
)
