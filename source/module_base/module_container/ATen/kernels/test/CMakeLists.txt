AddTest(
  TARGET container_kernels_uts
  LIBS ${math_libs}
  SOURCES blas_test.cpp lapack_test.cpp
  memory_test.cpp linalg_test.cpp
)

target_link_libraries(container_kernels_uts container base device)