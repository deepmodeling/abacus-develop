# list(APPEND module_device_srcs
#   memory_op.cpp
#   device.cpp
# )

# if(USE_CUDA)
#   list(APPEND module_device_srcs
#     cuda/memory_op.cu
#   )
# endif()


# if(USE_ROCM)
#   hip_add_library(module_device_rocm STATIC
#     rocm/memory_op.hip.cu
#   )
#   target_link_libraries(
#     device 
#     module_device_rocm 
#     hip::host
#     hip::device
#     hip::hipfft
#     roc::hipblas
#     roc::hipsolver
#   )
# endif()

# add_library(device OBJECT ${device_srcs})

# if(USE_CUDA)
#   target_link_libraries(
#     device 
#   )
# elseif(USE_ROCM)
#   target_link_libraries(
#     device 
#     device_rocm 
#     hip::host
#     hip::device
#     hip::hipfft
#     roc::hipblas
#     roc::hipsolver
#   )
# endif()