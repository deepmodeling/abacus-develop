add_library(
  hsolver_linalg OBJECT
  orb_buffer.cpp
  exp_action.cpp
  exp_taylor.cpp
  exp_lanczos.cpp
  hsolver_linear.cpp
  linear_bicgstab.cpp
  linear_cgs.cpp
  linear_gmres.cpp
  linear_workspace.cpp
  linear_algebra.cpp
  linear_low_rank.cpp
  cn_subspace.cpp
  orthonormal.cpp
  orthonormal_small.cpp)

if(ENABLE_COVERAGE)
  add_coverage(hsolver_linalg)
endif()

if(BUILD_TESTING)
  add_subdirectory(unittests)
endif()
