set(_adapter_test_main)
if(ENABLE_MPI)
  list(APPEND _adapter_test_main ../../source_base/test/mpi_test_main.cpp)
endif()

addtest(
  TARGET MODULE_CELL_mdcell_adapter
  LIBS parameter base device cell domain_decomposition neighbor symmetry
  SOURCES test_mdcell_adapter.cpp
          ${_adapter_test_main}
          ../module_neighlist/neighbor_search.cpp
          ../module_neighlist/bin_manager.cpp
          ../module_neighlist/page_allocator.cpp)

if(ENABLE_MPI)
  add_test(
    NAME MODULE_CELL_mdcell_adapter_np2
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
            $<TARGET_FILE:MODULE_CELL_mdcell_adapter>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(MODULE_CELL_mdcell_adapter_np2
                       PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" TIMEOUT 60)
  add_test(
    NAME MODULE_CELL_mdcell_adapter_np4
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
            $<TARGET_FILE:MODULE_CELL_mdcell_adapter>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(MODULE_CELL_mdcell_adapter_np4
                       PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" TIMEOUT 60)
endif()

addtest(
  TARGET MODULE_CELL_index_test
  LIBS base device
  SOURCES test_cell_index.cpp ../cell_index.cpp)

addtest(
  TARGET MODULE_CELL_klist_validation_test
  LIBS parameter base device cell domain_decomposition neighbor symmetry
  SOURCES test_klist.cpp
          ../module_neighlist/neighbor_search.cpp
          ../module_neighlist/bin_manager.cpp
          ../module_neighlist/page_allocator.cpp)

addtest(
  TARGET MODULE_CELL_read_species_malformed_test
  LIBS parameter base device cell domain_decomposition neighbor symmetry
  SOURCES test_read_atom_species.cpp
          ../module_neighlist/neighbor_search.cpp
          ../module_neighlist/bin_manager.cpp
          ../module_neighlist/page_allocator.cpp)

addtest(
  TARGET MODULE_CELL_read_atoms_malformed_test
  LIBS parameter base device cell domain_decomposition neighbor symmetry
  SOURCES test_read_atoms.cpp
          ../module_neighlist/neighbor_search.cpp
          ../module_neighlist/bin_manager.cpp
          ../module_neighlist/page_allocator.cpp)
