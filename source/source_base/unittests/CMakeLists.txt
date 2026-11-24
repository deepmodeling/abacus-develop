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

# The sub-communicator helpers are MPI-only, so this target keeps __MPI even
# though the directory disables it for the others, and the test is registered
# only where real collectives are available.
if(ENABLE_MPI)
  addtest(
    TARGET MODULE_BASE_parallel_common
    LIBS parameter base device
    KEEP_FEATURE_DEFINITIONS __MPI
    SOURCES test_parallel_common.cpp)
  add_test(
    NAME MODULE_BASE_parallel_common_4
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:MODULE_BASE_parallel_common> ${MPIEXEC_POSTFLAGS}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(MODULE_BASE_parallel_common_4 PROPERTIES PROCESSORS 4
                                                                 TIMEOUT 60)
endif()
