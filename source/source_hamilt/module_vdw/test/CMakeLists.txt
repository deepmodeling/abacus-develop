abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)

install(FILES c6.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES r0.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES r0_zero.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

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
