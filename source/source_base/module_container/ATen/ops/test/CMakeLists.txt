AddTest(
  TARGET MODULE_BASE_container_ops_uts
  LIBS parameter ${math_libs}
  SOURCES einsum_op_test.cpp linalg_op_test.cpp ../../kernels/lapack.cpp
)

target_link_libraries(MODULE_BASE_container_ops_uts container base device)
