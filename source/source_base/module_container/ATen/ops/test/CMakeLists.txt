AddTest(
  TARGET MODULE_BASE_container_ops_uts
  LIBS parameter
  SOURCES einsum_op_test.cpp linalg_op_test.cpp ../../kernels/lapack.cpp
)

target_link_libraries(MODULE_BASE_container_ops_uts PRIVATE container base device)
