#add_subdirectory(kernels/cuda)
if(ENABLE_LCAO)

list(APPEND objects
    gint.cpp
    gint_gamma_env.cpp
    gint_gamma_vl.cpp
    gint_fvl.cpp
    gint_rho.cpp
    gint_tau.cpp
    gint_vl.cpp
    gint_k_env.cpp
    gint_k_sparse1.cpp
    gint_k_pvpr.cpp
    gint_k_pvdpr.cpp
    gint_tools.cpp
    grid_bigcell.cpp
    grid_meshball.cpp
    grid_meshcell.cpp
    grid_meshk.cpp
    grid_technique.cpp
    gint_force_cpu_interface.cpp
    gint_rho_cpu_interface.cpp
    gint_vl_cpu_interface.cpp
    cal_psir_ylm.cpp
    cal_dpsir_ylm.cpp
    cal_ddpsir_ylm.cpp
    mult_psi_dmr.cpp
    init_orb.cpp
)

if(USE_CUDA)
  list(APPEND objects
      gint_gpu_interface.cpp
      kernels/cuda/cuda_tools.cu
      kernels/cuda/gint_vl.cu
      kernels/cuda/gint_rho.cu
      kernels/cuda/gint_force.cu
      gint_vl_gpu.cu
      gint_rho_gpu.cu
      gint_force_gpu.cu
      kernels/cuda/gemm_selector.cu
      kernels/cuda/code_gen_00.cu
      kernels/cuda/code_gen_01.cu
      kernels/cuda/code_gen_02.cu
      kernels/cuda/code_gen_03.cu
      kernels/cuda/code_gen_04.cu
      kernels/cuda/code_gen_05.cu
      kernels/cuda/code_gen_06.cu
      kernels/cuda/code_gen_07.cu
      kernels/cuda/code_gen_08.cu
      kernels/cuda/code_gen_09.cu
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

IF (BUILD_TESTING)
  if(ENABLE_MPI)
    add_subdirectory(test)
  endif()
endif()

endif()