AddTest(
  TARGET test_grad_matrix_degenerate
  LIBS base parameter ${math_libs}
  SOURCES test_grad_matrix_degenerate.cpp ../grad_matrix_degenerate.cpp
)
