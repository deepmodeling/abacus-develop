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
      kernels/cuda/gint_vl.cu
      kernels/cuda/vbatch_matrix_multiple/cuda_tools.cu
      kernels/cuda/gtask_vl.cpp
      kernels/cuda/gint_vl_gpu.cu
      kernels/cuda/vbatch_matrix_multiple/vbatch_matrix_mul.cu
      kernels/cuda/gint_rho.cu
      kernels/cuda/gtask_rho.cpp
      kernels/cuda/gint_rho_gpu.cu
      kernels/cuda/gint_force_gpu.cu
      kernels/cuda/gint_force.cu
      kernels/cuda/gtask_force.cpp
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
