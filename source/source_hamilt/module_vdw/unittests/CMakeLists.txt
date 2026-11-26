abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)

abacus_copy_test_files(support)

addtest(
  TARGET MODULE_HAMILT_vdwTest
  LIBS parameter base device vdw
  SOURCES test_vdw.cpp)

if(ENABLE_DFTD4 AND ENABLE_MPI)
  add_test(
    NAME MODULE_HAMILT_vdwd4Test_2np
    COMMAND
      ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:MODULE_HAMILT_vdwTest> --gtest_filter=vdwd4Test.*
      ${MPIEXEC_POSTFLAGS}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif()
