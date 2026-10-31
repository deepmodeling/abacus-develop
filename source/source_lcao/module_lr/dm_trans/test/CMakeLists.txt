if(ENABLE_MPI)
  set(test_dm_diff_main ../../../../source_base/test/mpi_test_main.cpp)
endif()
abacus_disable_feature_definitions(__LIBXC)
addtest(
  TARGET MODULE_LR_dm_trans_test
  LIBS parameter psi base device container
  SOURCES dm_trans_test.cpp
          ../../utils/lr_util.cpp
          ../dm_trans_parallel.cpp
          ../dm_trans_serial.cpp
          # ../../../source_base/module_container/ATen/core/tensor.cpp
          # ../../../source_base/module_container/ATen/core/tensor_shape.cpp
          # ../../../source_base/module_container/ATen/core/tensor_buffer.cpp
          # ../../../source_base/module_container/ATen/core/tensor_map.cpp
          # ../../../source_base/module_container/ATen/core/tensor_types.cpp
          # ../../../source_base/module_container/base/core/cpu_allocator.cpp
          # ../../../source_base/module_container/base/core/refcount.cpp
          # ../../../source_base/module_container/ATen/kernels/memory_impl.cpp
)

addtest(
  TARGET MODULE_LR_dm_diff_test
  LIBS psi base parameter ${math_libs} device container
  SOURCES test_dm_diff.cpp ${test_dm_diff_main} ../dm_diff_ser.cpp
          ../dm_diff_par.cpp ../../utils/lr_util.cpp)
