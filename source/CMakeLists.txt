add_subdirectory(source_base)
add_subdirectory(source_cell)
add_subdirectory(source_psi)
add_subdirectory(source_estate)
add_subdirectory(source_hamilt)
add_subdirectory(source_pw)
add_subdirectory(source_lcao)
add_subdirectory(source_hsolver)
add_subdirectory(source_basis/module_ao)
add_subdirectory(source_basis/module_nao)
add_subdirectory(source_md)
add_subdirectory(source_basis/module_pw)
add_subdirectory(source_esolver)
add_subdirectory(source_lcao/module_gint)
add_subdirectory(source_io)
add_subdirectory(source_relax)
add_subdirectory(source_lcao/module_ri)
add_subdirectory(module_parameter)
add_subdirectory(module_lr)

# add by jghan
add_subdirectory(source_lcao/module_rdmft)

add_library(
    driver
    OBJECT
    source_main/driver.cpp
    source_main/driver_run.cpp
)

list(APPEND device_srcs
  source_pw/module_pwdft/kernels/nonlocal_op.cpp
  source_pw/module_pwdft/kernels/veff_op.cpp
  source_pw/module_pwdft/kernels/ekinetic_op.cpp
  source_pw/module_pwdft/kernels/meta_op.cpp
  source_pw/module_stodft/kernels/hpsi_norm_op.cpp
  source_basis/module_pw/kernels/pw_op.cpp
  source_hsolver/kernels/dngvd_op.cpp
  source_hsolver/kernels/bpcg_kernel_op.cpp
  source_estate/kernels/elecstate_op.cpp

  # source_psi/kernels/psi_memory_op.cpp
  # source_psi/kernels/device.cpp

  source_base/module_device/device.cpp
  source_base/module_device/memory_op.cpp
  source_base/kernels/math_kernel_op.cpp
  source_base/kernels/math_kernel_op_vec.cpp

  source_pw/module_pwdft/kernels/force_op.cpp
  source_pw/module_pwdft/kernels/stress_op.cpp
  source_pw/module_pwdft/kernels/onsite_op.cpp
  source_pw/module_pwdft/kernels/wf_op.cpp
  source_pw/module_pwdft/kernels/vnl_op.cpp
  source_base/kernels/math_ylm_op.cpp
  source_hamilt/module_xc/kernels/xc_functional_op.cpp
)

if(USE_CUDA)
  list(APPEND device_srcs
    source_pw/module_pwdft/kernels/cuda/nonlocal_op.cu
    source_pw/module_pwdft/kernels/cuda/veff_op.cu
    source_pw/module_pwdft/kernels/cuda/ekinetic_op.cu
    source_pw/module_pwdft/kernels/cuda/meta_op.cu
    source_pw/module_stodft/kernels/cuda/hpsi_norm_op.cu
    source_pw/module_pwdft/kernels/cuda/onsite_op.cu
    source_basis/module_pw/kernels/cuda/pw_op.cu
    source_hsolver/kernels/cuda/dngvd_op.cu
    source_hsolver/kernels/cuda/bpcg_kernel_op.cu
    source_estate/kernels/cuda/elecstate_op.cu

    # source_psi/kernels/cuda/memory_op.cu
    source_base/module_device/cuda/memory_op.cu

    source_pw/module_pwdft/kernels/cuda/force_op.cu
    source_pw/module_pwdft/kernels/cuda/stress_op.cu
    source_pw/module_pwdft/kernels/cuda/wf_op.cu
    source_pw/module_pwdft/kernels/cuda/vnl_op.cu
    source_base/kernels/cuda/math_ylm_op.cu
    source_base/kernels/cuda/math_kernel_op.cu
    source_base/kernels/cuda/math_kernel_op_vec.cu
    source_hamilt/module_xc/kernels/cuda/xc_functional_op.cu
  )
endif()

if(USE_ROCM)
  hip_add_library(device_rocm STATIC
    source_pw/module_pwdft/kernels/rocm/nonlocal_op.hip.cu
    source_pw/module_pwdft/kernels/rocm/veff_op.hip.cu
    source_pw/module_pwdft/kernels/rocm/ekinetic_op.hip.cu
    source_pw/module_pwdft/kernels/rocm/meta_op.hip.cu
    source_pw/module_pwdft/kernels/rocm/onsite_op.hip.cu
    source_pw/module_stodft/kernels/rocm/hpsi_norm_op.hip.cu
    source_basis/module_pw/kernels/rocm/pw_op.hip.cu
    source_hsolver/kernels/rocm/dngvd_op.hip.cu
    source_hsolver/kernels/rocm/bpcg_kernel_op.hip.cu
    source_estate/kernels/rocm/elecstate_op.hip.cu

    # source_psi/kernels/rocm/memory_op.hip.cu
    source_base/module_device/rocm/memory_op.hip.cu

    source_pw/module_pwdft/kernels/rocm/force_op.hip.cu
    source_pw/module_pwdft/kernels/rocm/stress_op.hip.cu
    source_pw/module_pwdft/kernels/rocm/wf_op.hip.cu
    source_pw/module_pwdft/kernels/rocm/vnl_op.hip.cu
    source_base/kernels/rocm/math_kernel_op.hip.cu
    source_base/kernels/rocm/math_kernel_op_vec.hip.cu
    source_base/kernels/rocm/math_ylm_op.hip.cu
    source_hamilt/module_xc/kernels/rocm/xc_functional_op.hip.cu
  )
endif()

if(USE_DSP)
  list(APPEND device_srcs
  source_base/kernels/dsp/dsp_connector.cpp
  )
endif()


add_library(device OBJECT ${device_srcs})

if(USE_CUDA)
  target_link_libraries(
    device 
    cusolver 
    cublas 
    cufft
  )
elseif(USE_ROCM)
  target_link_libraries(
    device 
    device_rocm 
    hip::host
    hip::device
    hip::hipfft
    roc::hipblas
    roc::hipsolver
  )
endif()

if(ENABLE_COVERAGE)
  add_coverage(driver)
endif()
