remove_definitions(-DUSE_LIBXC)
AddTest(
  TARGET dm_trans_test
  LIBS parameter psi base ${math_libs} device container
  SOURCES dm_trans_test.cpp ../../utils/lr_util.cpp ../dm_trans_parallel.cpp ../dm_trans_serial.cpp
  # ../../../source_base/module_container/ATen/core/tensor.cpp
  # ../../../source_base/module_container/ATen/core/tensor_shape.cpp
  # ../../../source_base/module_container/ATen/core/tensor_buffer.cpp
  # ../../../source_base/module_container/ATen/core/tensor_map.cpp
  # ../../../source_base/module_container/ATen/core/tensor_types.cpp
  # ../../../source_base/module_container/base/core/cpu_allocator.cpp
  # ../../../source_base/module_container/base/core/refcount.cpp
  # ../../../source_base/module_container/ATen/kernels/memory_impl.cpp
)