abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)
abacus_disable_feature_definitions(__EXX)

set(TD_LINEAR_TEST_SOURCES
  ../hsolver_linear.cpp
  ../linear_bicgstab.cpp
  ../linear_cgs.cpp
  ../linear_gmres.cpp
  ../linear_algebra.cpp
  ../cn_subspace.cpp
  ../linear_low_rank.cpp
  ../linear_workspace.cpp
)
AddTest(
  TARGET MODULE_HSOLVER_linear
  LIBS parameter base device container
  SOURCES test_hsolver_linear.cpp test_cn_subspace.cpp test_linear_low_rank.cpp ${TD_LINEAR_TEST_SOURCES}
)
AddTest(
  TARGET MODULE_HSOLVER_pw_tddft
  LIBS parameter base device container psi planewave
  SOURCES test_cn_solver_pw.cpp test_orthonormal.cpp test_exp_action.cpp
    ../exp_action.cpp ../exp_lanczos.cpp ../exp_taylor.cpp
    ../cn_solver_pw.cpp ../td_orbital_ops.cpp ../orb_buffer.cpp ../orthonormal.cpp ../orthonormal_small.cpp ${TD_LINEAR_TEST_SOURCES}
)
if(ENABLE_MPI)
  add_test(
    NAME MODULE_HSOLVER_linear_parallel
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:MODULE_HSOLVER_linear> ${MPIEXEC_POSTFLAGS}
  )
  set_tests_properties(MODULE_HSOLVER_linear_parallel PROPERTIES PROCESSORS 2 TIMEOUT 120)
endif()
