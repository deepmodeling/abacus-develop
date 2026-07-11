abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)
abacus_disable_feature_definitions(__EXX)



# Split tests: create one test target per test source to keep failures isolated.

AddTest(
  TARGET MODULE_CELL_NEIGHBOR_neighbor_search
  LIBS parameter base device
  SOURCES
    neighbor_search_test.cpp
  ../neighbor_search.cpp
  ../bin_manager.cpp
  ../page_allocator.cpp
  ../unitcell_lite.cpp
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
  add_executable(MODULE_CELL_NEIGHBOR_neighbor_search_mpi_benchmark
    neighbor_search_mpi_benchmark.cpp
    ../domain_decomposition.cpp
    ../neighbor_search.cpp
    ../bin_manager.cpp
    ../page_allocator.cpp
    ../unitcell_lite.cpp
    ../../../source_base/matrix.cpp
    ../../../source_base/matrix3.cpp
    ../../../source_base/tool_quit.cpp
  )
  target_include_directories(MODULE_CELL_NEIGHBOR_neighbor_search_mpi_benchmark PRIVATE ${ABACUS_SOURCE_DIR})
  target_compile_definitions(MODULE_CELL_NEIGHBOR_neighbor_search_mpi_benchmark PRIVATE __NORMAL)
  target_link_libraries(MODULE_CELL_NEIGHBOR_neighbor_search_mpi_benchmark
    PRIVATE
    Threads::Threads MPI::MPI_CXX
  )
  if(USE_OPENMP)
    target_link_libraries(MODULE_CELL_NEIGHBOR_neighbor_search_mpi_benchmark PRIVATE OpenMP::OpenMP_CXX)
  endif()
  install(TARGETS MODULE_CELL_NEIGHBOR_neighbor_search_mpi_benchmark DESTINATION ${CMAKE_BINARY_DIR}/tests)
  add_test(NAME MODULE_CELL_NEIGHBOR_neighbor_search_mpi_benchmark_np4
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
            $<TARGET_FILE:MODULE_CELL_NEIGHBOR_neighbor_search_mpi_benchmark>
            12 12 12 2 1.75 1.0 0.2 1
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
endif()
