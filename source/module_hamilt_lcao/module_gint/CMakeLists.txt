#add_subdirectory(kernels/cuda)

list(APPEND objects
    gint.cpp
    gint_gamma.cpp
    gint_gamma_env.cpp
    gint_gamma_vl.cpp
    gint_fvl.cpp
    gint_rho.cpp
    gint_tau.cpp
    gint_vl.cpp
    gint_k_env.cpp
    gint_k_sparse.cpp
    gint_k_sparse1.cpp
    gint_k_pvpr.cpp
    gint_k_pvdpr.cpp
    gint_tools.cpp
    grid_bigcell.cpp
    grid_meshball.cpp
    grid_meshcell.cpp
    grid_meshk.cpp
    grid_technique.cpp
)

if(USE_CUDA)
  list(APPEND objects
      kernels/cuda/vbatch_matrix_multiple/cuda_tools.cu
      kernels/cuda/vbatch_matrix_multiple/vbatch_matrix_mul.cu
      kernels/cuda/gint_vl.cu
      kernels/cuda/gint_rho.cu
      kernels/cuda/gint_force.cu
      gint_vl_gpu.cu
      gint_rho_gpu.cu
      gint_force_gpu.cu
      gtask_vl.cpp
      gtask_rho.cpp
      gtask_force.cpp
  )
endif()

add_library(
    gint
    OBJECT
    ${objects}
)

if(ENABLE_COVERAGE)
  add_coverage(gint)
endif()
