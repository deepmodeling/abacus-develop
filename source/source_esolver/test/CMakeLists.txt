remove_definitions(-D__MPI)
remove_definitions(-D__LCAO)

install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

AddTest(
  TARGET MODULE_ESOLVER_esolver_dp_test
  LIBS parameter ${math_libs} base device 
  SOURCES esolver_dp_test.cpp ../esolver_dp.cpp ../../source_io/module_output/cif_io.cpp ../../source_io/module_output/output_log.cpp
)
