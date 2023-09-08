add_subdirectory(module_base)
add_subdirectory(module_cell)
add_subdirectory(module_psi)
add_subdirectory(module_elecstate)
add_subdirectory(module_hamilt_general)
add_subdirectory(module_hamilt_pw)
add_subdirectory(module_hamilt_lcao)
add_subdirectory(module_hsolver)
add_subdirectory(module_basis/module_ao)
add_subdirectory(module_basis/module_nao)
add_subdirectory(module_md)
add_subdirectory(module_basis/module_pw)
add_subdirectory(module_esolver)
add_subdirectory(module_hamilt_lcao/module_gint)
add_subdirectory(module_io)
add_subdirectory(module_relax)
add_subdirectory(module_ri)

add_library(
    driver
    OBJECT
    driver.cpp
    driver_run.cpp
)

list(APPEND device_srcs
  module_hamilt_pw/hamilt_pwdft/kernels/nonlocal_op.cpp
  module_hamilt_pw/hamilt_pwdft/kernels/veff_op.cpp
  module_hamilt_pw/hamilt_pwdft/kernels/ekinetic_op.cpp
  module_hamilt_pw/hamilt_pwdft/kernels/meta_op.cpp
  module_basis/module_pw/kernels/pw_op.cpp
  module_hsolver/kernels/dngvd_op.cpp
  module_hsolver/kernels/math_kernel_op.cpp
  module_elecstate/kernels/elecstate_op.cpp
  module_psi/kernels/psi_memory_op.cpp
  module_psi/kernels/device.cpp
  module_hamilt_pw/hamilt_pwdft/kernels/force_op.cpp
  module_hamilt_pw/hamilt_pwdft/kernels/stress_op.cpp
  module_hamilt_pw/hamilt_pwdft/kernels/wf_op.cpp
  module_hamilt_pw/hamilt_pwdft/kernels/vnl_op.cpp
  module_base/kernels/math_op.cpp
)

if(USE_CUDA)
  list(APPEND device_srcs
    module_hamilt_pw/hamilt_pwdft/kernels/cuda/nonlocal_op.cu
    module_hamilt_pw/hamilt_pwdft/kernels/cuda/veff_op.cu
    module_hamilt_pw/hamilt_pwdft/kernels/cuda/ekinetic_op.cu
    module_hamilt_pw/hamilt_pwdft/kernels/cuda/meta_op.cu
    module_basis/module_pw/kernels/cuda/pw_op.cu
    module_hsolver/kernels/cuda/dngvd_op.cu
    module_hsolver/kernels/cuda/math_kernel_op.cu
    module_elecstate/kernels/cuda/elecstate_op.cu
    module_psi/kernels/cuda/memory_op.cu
    module_hamilt_pw/hamilt_pwdft/kernels/cuda/force_op.cu
    module_hamilt_pw/hamilt_pwdft/kernels/cuda/stress_op.cu
    module_hamilt_pw/hamilt_pwdft/kernels/cuda/wf_op.cu
    module_hamilt_pw/hamilt_pwdft/kernels/cuda/vnl_op.cu
    module_base/kernels/cuda/math_op.cu
  )
endif()

if(USE_ROCM)
  hip_add_library(device_rocm STATIC
    module_hamilt_pw/hamilt_pwdft/kernels/rocm/nonlocal_op.hip.cu
    module_hamilt_pw/hamilt_pwdft/kernels/rocm/veff_op.hip.cu
    module_hamilt_pw/hamilt_pwdft/kernels/rocm/ekinetic_op.hip.cu
    module_hamilt_pw/hamilt_pwdft/kernels/rocm/meta_op.hip.cu
    module_basis/module_pw/kernels/rocm/pw_op.hip.cu
    module_hsolver/kernels/rocm/dngvd_op.hip.cu
    module_hsolver/kernels/rocm/math_kernel_op.hip.cu
    module_elecstate/kernels/rocm/elecstate_op.hip.cu
    module_psi/kernels/rocm/memory_op.hip.cu
    module_hamilt_pw/hamilt_pwdft/kernels/rocm/force_op.hip.cu
    module_hamilt_pw/hamilt_pwdft/kernels/rocm/stress_op.hip.cu
    module_hamilt_pw/hamilt_pwdft/kernels/rocm/wf_op.hip.cu
    module_hamilt_pw/hamilt_pwdft/kernels/rocm/vnl_op.hip.cu
    module_base/kernels/rocm/math_op.hip.cu
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
  set_property(
    TARGET device
    PROPERTY CUDA_ARCHITECTURES
    60 # P100
    70 # V100
    75 # T4
    80 # A100
  )
elseif(USE_ROCM)
  target_link_libraries(
    device 
    device_rocm 
    hip::device
    hip::host
    hip::hipfft
    roc::hipblas
  )
endif()

if(ENABLE_COVERAGE)
  add_coverage(driver)
endif()
