abacus_disable_feature_definitions(__MPI)
abacus_copy_test_files(support)

addtest(
  TARGET MODULE_BASE_csr_reader
  LIBS parameter base device
  SOURCES test_csr_file_reader.cpp test_csr_reader.cpp)

addtest(
  TARGET MODULE_BASE_filename
  LIBS parameter base device
  SOURCES test_filename.cpp)

addtest(
  TARGET MODULE_BASE_integral
  LIBS parameter base device
  SOURCES test_math_integral.cpp)

addtest(
  TARGET MODULE_BASE_net_force
  LIBS parameter base device
  SOURCES test_mathzone.cpp)

# alltoallv_data is an MPI wrapper, so this target keeps __MPI.
addtest(
  TARGET MODULE_BASE_test_parallel_device
  LIBS MPI::MPI_CXX base device parameter
  SOURCES test_parallel_device.cpp ../test/mpi_test_main.cpp
  KEEP_FEATURE_DEFINITIONS __MPI)

add_test(
  NAME MODULE_BASE_test_parallel_device_np4
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:MODULE_BASE_test_parallel_device> ${MPIEXEC_POSTFLAGS})
set_tests_properties(MODULE_BASE_test_parallel_device_np4
                     PROPERTIES PROCESSORS 4 TIMEOUT 120)
