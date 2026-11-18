addtest(
  TARGET MODULE_CELL_mdcell_adapter
  LIBS parameter base device cell domain_decomposition neighbor symmetry
  SOURCES test_mdcell_adapter.cpp
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
endif()
