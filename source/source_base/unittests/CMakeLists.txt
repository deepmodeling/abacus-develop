abacus_disable_feature_definitions(__MPI)
abacus_copy_test_files(support)

addtest(
  TARGET MODULE_BASE_csr_reader
  LIBS parameter base device
  SOURCES test_csr_file_reader.cpp test_csr_reader.cpp)

AddTest(
  TARGET MODULE_BASE_filename
  LIBS parameter base device
  SOURCES test_filename.cpp
)

# The sub-communicator helpers need real collectives, so this target keeps
# __MPI even though the directory disables it for the others. It is built
# directly rather than through addtest(): that helper also registers a
# single-process ctest entry, and this test requires four ranks.
if(ENABLE_MPI)
  add_executable(MODULE_BASE_parallel_common test_parallel_common.cpp)
  # The test drives a plain `main` around MPI_Init and its own require() checks
  # rather than a gtest fixture, so it needs no gtest/gmock linkage.
  target_link_libraries(MODULE_BASE_parallel_common
                        PRIVATE parameter base device abacus::linalg_libs)
  target_compile_definitions(MODULE_BASE_parallel_common PRIVATE __MPI)
  add_test(
    NAME MODULE_BASE_parallel_common_4
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:MODULE_BASE_parallel_common> ${MPIEXEC_POSTFLAGS}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  # No per-case time limit: ctest's own default applies, so a slow or
  # oversubscribed CI machine cannot turn a passing test into a timeout.
  set_tests_properties(MODULE_BASE_parallel_common_4 PROPERTIES PROCESSORS 4)
endif()
