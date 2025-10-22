#add_subdirectory(kernels/cuda)
if(ENABLE_LCAO)

list(APPEND objects
    biggrid_info.cpp
    big_grid.cpp
    divide_info.cpp
    gint_atom.cpp
    gint_info.cpp
    gint.cpp
    gint_vl.cpp
    gint_vl_metagga.cpp
    gint_vl_nspin4.cpp
    gint_vl_metagga_nspin4.cpp
    gint_rho.cpp
    gint_tau.cpp
    gint_fvl.cpp
    gint_fvl_meta.cpp
    gint_env_gamma.cpp
    gint_env_k.cpp
    gint_dvlocal.cpp
    localcell_info.cpp
    phi_operator.cpp
    set_ddphi.cpp
    unitcell_info.cpp
    gint_common.cpp
    gint_interface.cpp
    )
if(USE_CUDA)
  list(APPEND objects
      kernel/gint_gpu_vars.cpp
      kernel/phi_operator_gpu.cu
      kernel/phi_operator_kernel.cu
      kernel/set_const_mem.cu
      batch_biggrid.cpp
      gint_vl_gpu.cpp
      gint_rho_gpu.cpp
      gint_fvl_gpu.cpp
      gint_vl_metagga_gpu.cpp
      gint_vl_nspin4_gpu.cpp
      gint_vl_metagga_nspin4_gpu.cpp
      gint_tau_gpu.cpp
      gint_fvl_meta_gpu.cpp
      kernel/dgemm_vbatch.cu
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

endif()