abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)

abacus_copy_test_files(
  c6.txt
  r0.txt
  r0_zero.txt)

AddTest(
  TARGET MODULE_HAMILT_vdwTest
  LIBS parameter base device vdw
  SOURCES vdw_test.cpp vdwd3_evaluator_test.cpp
)

if(ENABLE_DFTD4 AND ENABLE_MPI)
  add_test(
    NAME MODULE_HAMILT_vdwd4Test_2np
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
            ${MPIEXEC_PREFLAGS} $<TARGET_FILE:MODULE_HAMILT_vdwTest>
            --gtest_filter=vdwd4Test.* ${MPIEXEC_POSTFLAGS}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
endif()

if(ENABLE_LIBVDWXC AND ENABLE_MPI)
  addtest(
    TARGET MODULE_HAMILT_nonlocal_xcTest
    LIBS xc_
         planewave
         parameter
         MPI::MPI_CXX
         Libxc::xc
         psi
         device
         container
         base
    SOURCES test_nonlocal_xc.cpp)
  set_tests_properties(MODULE_HAMILT_nonlocal_xcTest
                       PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")
  add_test(
    NAME MODULE_HAMILT_nonlocal_xcTest_2np
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:MODULE_HAMILT_nonlocal_xcTest> ${MPIEXEC_POSTFLAGS}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(MODULE_HAMILT_nonlocal_xcTest_2np
                       PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" PROCESSORS 2)
  # four processes give an uneven z distribution and an empty libvdwxc slab
  add_test(
    NAME MODULE_HAMILT_nonlocal_xcTest_4np
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:MODULE_HAMILT_nonlocal_xcTest> ${MPIEXEC_POSTFLAGS}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(MODULE_HAMILT_nonlocal_xcTest_4np
                       PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" PROCESSORS 4)
endif()
