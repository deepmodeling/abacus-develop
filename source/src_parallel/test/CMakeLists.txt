AddTest(
  TARGET ParaCommon
  LIBS MPI::MPI_CXX 
  SOURCES parallel_common_test.cpp ../../module_base/global_variable.cpp ../parallel_common.cpp
)

AddTest(
  TARGET ParaGlobal
  LIBS MPI::MPI_CXX 
  SOURCES parallel_global_test.cpp ../../module_base/global_variable.cpp ../parallel_global.cpp
)

AddTest(
  TARGET ParaKpoints
  LIBS MPI::MPI_CXX
  SOURCES parallel_kpoints_test.cpp ../../module_base/global_variable.cpp ../parallel_global.cpp ../parallel_common.cpp ../parallel_kpoints.cpp
)

install(FILES parallel_common_test.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES parallel_global_test.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES parallel_kpoints_test.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

find_program(BASH bash)

add_test(NAME parallel_common_test
      COMMAND ${BASH} parallel_common_test.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME parallel_global_test
      COMMAND ${BASH} parallel_global_test.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME parallel_kpoints_test
      COMMAND ${BASH} parallel_kpoints_test.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
