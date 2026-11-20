addtest(
  TARGET MODULE_BASE_test_parallel_device
  LIBS MPI::MPI_CXX base device parameter
  SOURCES test_parallel_device.cpp ../test/mpi_test_main.cpp)

add_test(
  NAME MODULE_BASE_test_parallel_device_np4
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:MODULE_BASE_test_parallel_device> ${MPIEXEC_POSTFLAGS}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(MODULE_BASE_test_parallel_device_np4
                     PROPERTIES PROCESSORS 4 TIMEOUT 120)
