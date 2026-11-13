abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)

addtest(
  TARGET MODULE_HAMILT_test_vdw
  LIBS parameter base device vdw
  SOURCES test_vdw.cpp)

if(ENABLE_LIBVDWXC AND ENABLE_MPI)
  set(nonlocal_xc_test_libs
      xc_
      planewave
      parameter
      MPI::MPI_CXX
      psi
      device
      container
      base)
  if(ENABLE_LIBXC)
    list(APPEND nonlocal_xc_test_libs Libxc::xc)
  endif()
  addtest(
    TARGET MODULE_HAMILT_nonlocal_xcTest
    LIBS ${nonlocal_xc_test_libs}
    SOURCES test_nonlocal_xc.cpp test_rvv10_workspace.cpp)
  # a collective-path regression shows as a hang, so every registration gets a timeout
  set_tests_properties(MODULE_HAMILT_nonlocal_xcTest
                       PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" TIMEOUT 300)
  add_test(
    NAME MODULE_HAMILT_nonlocal_xcTest_2np
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:MODULE_HAMILT_nonlocal_xcTest> ${MPIEXEC_POSTFLAGS}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(
    MODULE_HAMILT_nonlocal_xcTest_2np
    PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" PROCESSORS 2 TIMEOUT 300)
  # the same with two threads: the element-wise loops of the adapter are OpenMP-parallel
  add_test(
    NAME MODULE_HAMILT_nonlocal_xcTest_2np_2threads
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:MODULE_HAMILT_nonlocal_xcTest> ${MPIEXEC_POSTFLAGS}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(
    MODULE_HAMILT_nonlocal_xcTest_2np_2threads
    PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2" PROCESSORS 4 TIMEOUT 300)
  # four processes give an uneven z distribution, an empty libvdwxc slab and, on the big-cell grid, a
  # process without any point
  add_test(
    NAME MODULE_HAMILT_nonlocal_xcTest_4np
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:MODULE_HAMILT_nonlocal_xcTest> ${MPIEXEC_POSTFLAGS}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(
    MODULE_HAMILT_nonlocal_xcTest_4np
    PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" PROCESSORS 4 TIMEOUT 300)
endif()
