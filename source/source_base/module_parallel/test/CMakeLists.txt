abacus_disable_feature_definitions(__MPI)
AddTest(
  TARGET MODULE_BASE_para_world
  SOURCES para_world_test.cpp ../para_world.cpp
)

AddTest(
  TARGET MODULE_BASE_para_collection
  SOURCES para_collection_test.cpp ../para_collection.cpp ../para_world.cpp ../para_kmesh_world.cpp
)

AddTest(
  TARGET MODULE_BASE_para_kmesh_world
  SOURCES para_kmesh_world_test.cpp ../para_kmesh_world.cpp ../para_world.cpp
)

AddTest(
  TARGET MODULE_BASE_para_pw_world
  SOURCES para_pw_world_test.cpp ../para_pw_world.cpp ../para_world.cpp
)

AddTest(
  TARGET MODULE_BASE_para_diag_world
  SOURCES para_diag_world_test.cpp ../para_diag_world.cpp ../para_world.cpp
)

AddTest(
  TARGET MODULE_BASE_para_rgrid_world
  SOURCES para_rgrid_world_test.cpp ../para_rgrid_world.cpp ../para_world.cpp
)

AddTest(
  TARGET MODULE_BASE_para_bdiff_ksame_world
  SOURCES para_bdiff_ksame_world_test.cpp ../para_bdiff_ksame_world.cpp ../para_world.cpp
)

AddTest(
  TARGET MODULE_BASE_para_matrix_world
  SOURCES para_matrix_world_test.cpp ../para_matrix_world.cpp ../para_world.cpp
)

AddTest(
  TARGET MODULE_BASE_para_mpi_func
  SOURCES para_mpi_func_test.cpp ../para_mpi_func.cpp ../para_world.cpp
)

AddTest(
  TARGET MODULE_BASE_para_setup
  SOURCES para_setup_test.cpp ../para_setup.cpp ../para_world.cpp ../para_collection.cpp
)

AddTest(
  TARGET MODULE_BASE_para_world_mpi
  LIBS MPI::MPI_CXX
  SOURCES para_world_mpi_test.cpp ../para_world.cpp
)
target_compile_definitions(MODULE_BASE_para_world_mpi PRIVATE __MPI)

AddTest(
  TARGET MODULE_BASE_para_collection_mpi
  LIBS MPI::MPI_CXX
  SOURCES para_collection_mpi_test.cpp ../para_collection.cpp ../para_world.cpp ../para_kmesh_world.cpp
)
target_compile_definitions(MODULE_BASE_para_collection_mpi PRIVATE __MPI)

AddTest(
  TARGET MODULE_BASE_para_mpi_func_mpi
  LIBS MPI::MPI_CXX
  SOURCES para_mpi_func_mpi_test.cpp ../para_mpi_func.cpp ../para_world.cpp
)
target_compile_definitions(MODULE_BASE_para_mpi_func_mpi PRIVATE __MPI)

AddTest(
  TARGET MODULE_BASE_para_setup_mpi
  LIBS MPI::MPI_CXX
  SOURCES para_setup_mpi_test.cpp ../para_setup.cpp ../para_world.cpp ../para_collection.cpp
)
target_compile_definitions(MODULE_BASE_para_setup_mpi PRIVATE __MPI)

# Built with add_executable (not AddTest) so that no direct-run CTest entry is
# created; the binary is only exercised through mpirun by the .sh test below,
# matching the multi-process requirement of these cases.
add_executable(MODULE_BASE_para_kmesh_world_mpi test_para_kmesh_world_mpi.cpp ../para_bdiff_ksame_world.cpp ../para_kmesh_world.cpp ../para_world.cpp)
target_link_libraries(MODULE_BASE_para_kmesh_world_mpi PRIVATE MPI::MPI_CXX GTest::gtest GTest::gtest_main abacus::linalg_libs)
target_compile_definitions(MODULE_BASE_para_kmesh_world_mpi PRIVATE __MPI)

file(COPY para_world_mpi_test.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY para_collection_mpi_test.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY para_mpi_func_mpi_test.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY para_setup_mpi_test.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY test_para_kmesh_world_mpi.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
find_program(BASH bash)
add_test(NAME MODULE_BASE_para_world_mpi_test
      COMMAND ${BASH} para_world_mpi_test.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME MODULE_BASE_para_collection_mpi_test
      COMMAND ${BASH} para_collection_mpi_test.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME MODULE_BASE_para_mpi_func_mpi_test
      COMMAND ${BASH} para_mpi_func_mpi_test.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME MODULE_BASE_para_setup_mpi_test
      COMMAND ${BASH} para_setup_mpi_test.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME MODULE_BASE_para_kmesh_world_mpi_test
      COMMAND ${BASH} test_para_kmesh_world_mpi.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
