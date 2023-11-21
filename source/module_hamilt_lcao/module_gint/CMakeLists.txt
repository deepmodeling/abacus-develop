#add_subdirectory(gint_gpu)

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
    gint_gpu/gint_vl.cu
    gint_gpu/vbatch_matrix_multiple/cuda_tools.cu
    gint_gpu/gpu_task_generator.cpp
    gint_gpu/gint_vl_gpu.cu
)

add_library(
    gint
    OBJECT
    ${objects}
)

if(ENABLE_COVERAGE)
  add_coverage(gint)
endif()
