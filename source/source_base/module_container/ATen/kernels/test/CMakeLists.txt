AddTest(
  TARGET MODULE_BASE_container_kernels_uts
  LIBS parameter ${math_libs}
  SOURCES blas_test.cpp lapack_test.cpp
  memory_test.cpp linalg_test.cpp
)

target_link_libraries(MODULE_BASE_container_kernels_uts container base device)
