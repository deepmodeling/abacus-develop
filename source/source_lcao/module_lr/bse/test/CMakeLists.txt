remove_definitions(-DUSE_LIBXC)
AddTest(
  TARGET MODULE_LR_bse_solver_test
  LIBS parameter base ${math_libs} ELPA::ELPA device psi 
  SOURCES hamilt_bse_solver_test.cpp ../hamilt_bse_solver.cpp
)