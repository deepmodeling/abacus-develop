add_subdirectory(module_base)
add_subdirectory(module_cell)
add_subdirectory(module_symmetry)
add_subdirectory(module_neighbor)
add_subdirectory(module_psi)
add_subdirectory(module_elecstate)
add_subdirectory(module_hamilt)
add_subdirectory(module_hsolver)
add_subdirectory(module_orbital)
add_subdirectory(module_md)
add_subdirectory(module_pw)
add_subdirectory(module_surchem)
add_subdirectory(module_deepks)
add_subdirectory(module_xc)
add_subdirectory(module_esolver)
add_subdirectory(module_gint)
add_subdirectory(src_io)
add_subdirectory(module_relax)
add_subdirectory(src_lcao)
add_subdirectory(src_parallel)
add_subdirectory(src_pdiag)
add_subdirectory(src_pw)
add_subdirectory(src_ri)
add_subdirectory(module_rpa)
add_subdirectory(module_vdw)
add_subdirectory(module_ri)

add_library(
    driver
    OBJECT
    driver.cpp
    driver_run.cpp
    input.cpp
    input_conv.cpp
)

list(APPEND device_srcs
    module_hamilt/src/nonlocal.cpp
    module_hamilt/src/veff.cpp
    module_hamilt/src/ekinetic.cpp
    module_pw/src/pw_multi_device.cpp
    module_hsolver/src/dngvd_op.cpp
    module_hsolver/src/math_kernel.cpp
    module_elecstate/src/elecstate_multi_device.cpp
    module_psi/src/memory_psi.cpp
    module_psi/src/device.cpp
    src_pw/src/force_multi_device.cpp
    src_pw/src/stress_multi_device.cpp
    module_base/src/math_multi_device.cpp
)

if(USE_CUDA)
  list(APPEND device_srcs
      module_hamilt/src/cuda/nonlocal.cu
      module_hamilt/src/cuda/veff.cu
      module_hamilt/src/cuda/ekinetic.cu
      module_pw/src/cuda/pw_multi_device.cu
      module_hsolver/src/cuda/dngvd_op.cu
      module_hsolver/src/cuda/math_kernel.cu
      module_elecstate/src/cuda/elecstate_multi_device.cu
      module_psi/src/cuda/memory.cu
      src_pw/src/cuda/force_multi_device.cu
      src_pw/src/cuda/stress_multi_device.cu
      module_base/src/cuda/math_multi_device.cu
  )
endif()

add_library(device OBJECT ${device_srcs})

if(USE_CUDA)
  target_link_libraries(device cusolver cublas cufft)
endif()

if(ENABLE_COVERAGE)
  add_coverage(driver)
endif()
