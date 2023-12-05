AddTest(
  TARGET base_ParaCommon
  LIBS MPI::MPI_CXX 
  SOURCES parallel_common_test.cpp ../global_variable.cpp ../parallel_common.cpp
)

AddTest(
  TARGET base_ParaGlobal
  LIBS MPI::MPI_CXX 
  SOURCES parallel_global_test.cpp ../global_variable.cpp ../para_json.cpp ../parallel_global.cpp
)

AddTest(
  TARGET base_ParaReduce
  LIBS MPI::MPI_CXX
  SOURCES parallel_reduce_test.cpp ../global_variable.cpp ../para_json.cpp ../parallel_global.cpp ../parallel_common.cpp ../parallel_reduce.cpp
)

install(FILES parallel_common_test.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES parallel_global_test.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES parallel_reduce_test.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

find_program(BASH bash)
add_test(NAME base_parallel_common_test
      COMMAND ${BASH} parallel_common_test.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME base_parallel_global_test
      COMMAND ${BASH} parallel_global_test.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME base_parallel_reduce_test
      COMMAND ${BASH} parallel_reduce_test.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
