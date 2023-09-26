AddTest(
  TARGET Container_Kernels_UTs
  LIBS ${math_libs}
  SOURCES blas_op_test.cpp lapack_op_test.cpp einsum_op_test.cpp
  linalg_op_test.cpp memory_op_test.cpp
)

target_link_libraries(Container_Kernels_UTs container base device)