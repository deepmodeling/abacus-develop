abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)
abacus_disable_feature_definitions(__EXX)



# Split tests: create one test target per test source to keep failures isolated.

AddTest(
  TARGET MODULE_CELL_NEIGHBOR_neighbor_search
  LIBS parameter base device cell symmetry
  SOURCES
    neighbor_search_test.cpp
  ../neighbor_search.cpp
  ../domain_decomposition.cpp
  ../bin_manager.cpp
  ../page_allocator.cpp
)

AddTest(
  TARGET MODULE_CELL_NEIGHBOR_bin_manager
  LIBS parameter base device
  SOURCES
    bin_manager_test.cpp
    ../bin_manager.cpp
    ../page_allocator.cpp
)

AddTest(
  TARGET MODULE_CELL_NEIGHBOR_allocator_and_list
  LIBS parameter base device
  SOURCES
  neighbor_list_test.cpp
  ../page_allocator.cpp
)

if(ENABLE_MPI)
  add_executable(MODULE_CELL_NEIGHBOR_mdcell_migrate_mpi
    md_cell_migrate_mpi_test.cpp
    ../../md_cell.cpp
    ../domain_decomposition.cpp
  )
  target_link_libraries(MODULE_CELL_NEIGHBOR_mdcell_migrate_mpi
    PRIVATE
    parameter base device MPI::MPI_CXX GTest::gtest_main GTest::gmock_main abacus::linalg_libs
  )
  install(TARGETS MODULE_CELL_NEIGHBOR_mdcell_migrate_mpi DESTINATION ${CMAKE_BINARY_DIR}/tests)
  add_test(NAME MODULE_CELL_NEIGHBOR_mdcell_migrate_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
            $<TARGET_FILE:MODULE_CELL_NEIGHBOR_mdcell_migrate_mpi>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

  add_executable(MODULE_CELL_NEIGHBOR_distributed_mdcell_reader
    distributed_mdcell_reader_test.cpp
    ../../distributed_mdcell_reader.cpp
    ../../md_cell.cpp
    ../domain_decomposition.cpp
    ../../../source_base/global_variable.cpp
    ../../../source_base/communication_domain.cpp
    ../../../source_base/matrix.cpp
    ../../../source_base/matrix3.cpp
    ../../../source_base/tool_quit.cpp
  )
  target_include_directories(MODULE_CELL_NEIGHBOR_distributed_mdcell_reader PRIVATE ${ABACUS_SOURCE_DIR})
  target_compile_definitions(MODULE_CELL_NEIGHBOR_distributed_mdcell_reader PRIVATE __NORMAL)
  target_link_libraries(MODULE_CELL_NEIGHBOR_distributed_mdcell_reader
    PRIVATE
    Threads::Threads MPI::MPI_CXX GTest::gtest GTest::gmock
  )
  install(TARGETS MODULE_CELL_NEIGHBOR_distributed_mdcell_reader DESTINATION ${CMAKE_BINARY_DIR}/tests)
  add_test(NAME MODULE_CELL_NEIGHBOR_distributed_mdcell_reader_np4
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
            $<TARGET_FILE:MODULE_CELL_NEIGHBOR_distributed_mdcell_reader>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

endif()
