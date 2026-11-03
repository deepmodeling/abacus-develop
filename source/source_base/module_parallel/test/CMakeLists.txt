abacus_disable_feature_definitions(__MPI)
addtest(TARGET MODULE_BASE_para_world SOURCES para_world_test.cpp
                                              ../para_world.cpp)

addtest(TARGET MODULE_BASE_para_collection
        SOURCES para_collection_test.cpp ../para_collection.cpp
                ../para_world.cpp ../para_kmesh_world.cpp)

addtest(TARGET MODULE_BASE_para_kmesh_world
        SOURCES para_kmesh_world_test.cpp ../para_kmesh_world.cpp
                ../para_world.cpp)

addtest(TARGET MODULE_BASE_para_pw_world
        SOURCES para_pw_world_test.cpp ../para_pw_world.cpp ../para_world.cpp)

addtest(TARGET MODULE_BASE_para_diag_world
        SOURCES para_diag_world_test.cpp ../para_diag_world.cpp
                ../para_world.cpp)

addtest(TARGET MODULE_BASE_para_rgrid_world
        SOURCES para_rgrid_world_test.cpp ../para_rgrid_world.cpp
                ../para_world.cpp)

addtest(TARGET MODULE_BASE_para_bgroup_world
        SOURCES para_bgroup_world_test.cpp ../para_bgroup_world.cpp
                ../para_world.cpp)

addtest(TARGET MODULE_BASE_para_matrix_world
        SOURCES para_matrix_world_test.cpp ../para_matrix_world.cpp
                ../para_world.cpp)

addtest(TARGET MODULE_BASE_para_mpi_func
        SOURCES para_mpi_func_test.cpp ../para_mpi_func.cpp ../para_world.cpp)

addtest(TARGET MODULE_BASE_para_setup
        SOURCES para_setup_test.cpp ../para_setup.cpp ../para_world.cpp
                ../para_collection.cpp)

addtest(
  TARGET MODULE_BASE_para_world_mpi
  LIBS MPI::MPI_CXX
  SOURCES para_world_mpi_test.cpp ../para_world.cpp)
target_compile_definitions(MODULE_BASE_para_world_mpi PRIVATE __MPI)

addtest(
  TARGET MODULE_BASE_para_collection_mpi
  LIBS MPI::MPI_CXX
  SOURCES para_collection_mpi_test.cpp ../para_collection.cpp ../para_world.cpp
          ../para_kmesh_world.cpp)
target_compile_definitions(MODULE_BASE_para_collection_mpi PRIVATE __MPI)

addtest(
  TARGET MODULE_BASE_para_mpi_func_mpi
  LIBS MPI::MPI_CXX
  SOURCES para_mpi_func_mpi_test.cpp ../para_mpi_func.cpp ../para_world.cpp)
target_compile_definitions(MODULE_BASE_para_mpi_func_mpi PRIVATE __MPI)

addtest(
  TARGET MODULE_BASE_para_setup_mpi
  LIBS MPI::MPI_CXX
  SOURCES para_setup_mpi_test.cpp ../para_setup.cpp ../para_world.cpp
          ../para_collection.cpp)
target_compile_definitions(MODULE_BASE_para_setup_mpi PRIVATE __MPI)

abacus_copy_test_files(para_world_mpi_test.sh para_collection_mpi_test.sh
                       para_mpi_func_mpi_test.sh para_setup_mpi_test.sh)
find_program(BASH bash)
add_test(
  NAME MODULE_BASE_para_world_mpi_test
  COMMAND ${BASH} para_world_mpi_test.sh
          "$<TARGET_FILE:MODULE_BASE_para_world_mpi>"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
add_test(
  NAME MODULE_BASE_para_collection_mpi_test
  COMMAND ${BASH} para_collection_mpi_test.sh
          "$<TARGET_FILE:MODULE_BASE_para_collection_mpi>"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
add_test(
  NAME MODULE_BASE_para_mpi_func_mpi_test
  COMMAND ${BASH} para_mpi_func_mpi_test.sh
          "$<TARGET_FILE:MODULE_BASE_para_mpi_func_mpi>"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
add_test(
  NAME MODULE_BASE_para_setup_mpi_test
  COMMAND ${BASH} para_setup_mpi_test.sh
          "$<TARGET_FILE:MODULE_BASE_para_setup_mpi>"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
