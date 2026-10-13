abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)

install(FILES c6.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES r0.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES r0_zero.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if(ENABLE_MPI)
  AddTest(
    TARGET MODULE_HAMILT_vdwTest
    LIBS parameter base device vdw
    SOURCES vdw_test.cpp vdwd3_evaluator_test.cpp
  )
endif()

if(ENABLE_DFTD4 AND ENABLE_MPI)
  add_test(
    NAME MODULE_HAMILT_vdwd4Test_2np
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
            ${MPIEXEC_PREFLAGS} $<TARGET_FILE:MODULE_HAMILT_vdwTest>
            --gtest_filter=vdwd4Test.* ${MPIEXEC_POSTFLAGS}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
endif()
if(ENABLE_LIBVDWXC AND NOT ENABLE_MPI AND NOT USE_CUDA AND NOT USE_ROCM)
  AddTest(
    TARGET MODULE_HAMILT_vdwTest_LIBVDWXC
    LIBS xc_ planewave base parameter device nonlocal_xc
    SOURCES test_nonlocal_xc_libvdwxc.cpp
  )
  if(ENABLE_LIBXC)
    target_link_libraries(MODULE_HAMILT_vdwTest_LIBVDWXC PRIVATE Libxc::xc)
  endif()
  target_link_libraries(MODULE_HAMILT_vdwTest_LIBVDWXC PRIVATE libvdwxc::vdwxc)
  set_tests_properties(MODULE_HAMILT_vdwTest_LIBVDWXC PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")
endif()

if(ENABLE_DFTD4 AND ENABLE_MPI)
  target_link_libraries(MODULE_HAMILT_vdwTest PRIVATE dftd4::dftd4)
endif()
