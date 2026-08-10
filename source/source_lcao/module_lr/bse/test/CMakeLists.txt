remove_definitions(-DUSE_LIBXC)
set(bse_solver_test_libs parameter base ${math_libs} device psi)
if(ENABLE_ELPA)
  list(APPEND bse_solver_test_libs ELPA::ELPA)
endif()

AddTest(
  TARGET MODULE_LR_bse_solver_test
  LIBS ${bse_solver_test_libs}
  SOURCES hamilt_bse_solver_test.cpp ../hamilt_bse_solver.cpp
)
