#add_subdirectory(kernels/cuda)
if(ENABLE_LCAO)

list(APPEND objects
    temp_gint/biggrid_info.cpp
    temp_gint/big_grid.cpp
    temp_gint/divide_info.cpp
    temp_gint/gint_atom.cpp
    temp_gint/gint_info.cpp
    temp_gint/gint.cpp
    temp_gint/gint_vl.cpp
    temp_gint/gint_vl_metagga.cpp
    temp_gint/gint_vl_nspin4.cpp
    temp_gint/gint_vl_metagga_nspin4.cpp
    temp_gint/gint_rho.cpp
    temp_gint/gint_tau.cpp
    temp_gint/gint_fvl.cpp
    temp_gint/gint_fvl_meta.cpp
    temp_gint/gint_env_gamma.cpp
    temp_gint/gint_env_k.cpp
    temp_gint/gint_dvlocal.cpp
    temp_gint/localcell_info.cpp
    temp_gint/phi_operator.cpp
    temp_gint/set_ddphi.cpp
    temp_gint/unitcell_info.cpp
    temp_gint/gint_common.cpp
    temp_gint/gint_interface.cpp
    )
if(USE_CUDA)
  list(APPEND objects
      temp_gint/kernel/gint_gpu_vars.cpp
      temp_gint/kernel/phi_operator_gpu.cu
      temp_gint/kernel/phi_operator_kernel.cu
      temp_gint/kernel/set_const_mem.cu
      temp_gint/batch_biggrid.cpp
      temp_gint/gint_vl_gpu.cpp
      temp_gint/gint_rho_gpu.cpp
      temp_gint/gint_fvl_gpu.cpp
      temp_gint/gint_vl_metagga_gpu.cpp
      temp_gint/gint_vl_nspin4_gpu.cpp
      temp_gint/gint_vl_metagga_nspin4_gpu.cpp
      temp_gint/gint_tau_gpu.cpp
      temp_gint/gint_fvl_meta_gpu.cpp
      temp_gint/kernel/dgemm_vbatch.cu
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