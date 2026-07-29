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
    gint_drho.cpp
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

if(USE_CUDA AND CUDAToolkit_VERSION VERSION_GREATER_EQUAL "12.8")
  # -static-global-template-stub exists since CUDA 12.8; its default flipped to
  # true in CUDA 13.0, giving the host-side stubs of explicitly instantiated
  # __global__ template kernels local linkage in whole-program mode (-rdc=false).
  # The instantiations in kernel/phi_operator_kernel.cu are called from
  # kernel/phi_operator_gpu.cu, so the link fails with "hidden symbol ... isn't
  # defined". Restore external linkage for the stubs instead of restructuring
  # the sources.
  target_compile_options(gint PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:-static-global-template-stub=false>)
endif()

if(BUILD_TESTING)
  add_subdirectory(test)
endif()

if(ENABLE_COVERAGE)
  add_coverage(gint)
endif()

endif()
