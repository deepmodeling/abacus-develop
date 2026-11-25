abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)
abacus_disable_feature_definitions(__EXX)

set(LINALG_TEST_SOURCES
    ../hsolver_linear.cpp
    ../linear_bicgstab.cpp
    ../linear_cgs.cpp
    ../linear_gmres.cpp
    ../linear_algebra.cpp
    ../cn_subspace.cpp
    ../linear_low_rank.cpp
    ../linear_workspace.cpp
    ../exp_action.cpp
    ../exp_lanczos.cpp
    ../exp_taylor.cpp
    ../orb_buffer.cpp
    ../orthonormal.cpp
    ../orthonormal_small.cpp)
addtest(
  TARGET MODULE_HSOLVER_linalg
  LIBS parameter base device container
  SOURCES test_hsolver_linear.cpp test_cn_subspace.cpp test_linear_low_rank.cpp
          test_orthonormal.cpp test_exp_action.cpp ${LINALG_TEST_SOURCES})
if(ENABLE_MPI)
  add_test(
    NAME MODULE_HSOLVER_linalg_parallel
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:MODULE_HSOLVER_linalg> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(MODULE_HSOLVER_linalg_parallel PROPERTIES PROCESSORS 2
                                                                 TIMEOUT 120)
endif()
