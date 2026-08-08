AddTest(
  TARGET MODULE_BASE_container_kernels_uts
  LIBS parameter
  SOURCES blas_test.cpp lapack_test.cpp
  memory_test.cpp linalg_test.cpp
)

target_link_libraries(MODULE_BASE_container_kernels_uts PRIVATE container base device)
