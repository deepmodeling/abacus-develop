// TODO: This is a temperary location for these functions.
// And will be moved to a global module(module base) later.
#ifndef MODULE_PSI_MEMORY_H_
#define MODULE_PSI_MEMORY_H_

#include <vector>
#include <stddef.h>
#include "module_psi/include/types.h"

namespace psi {
namespace memory {

template <typename FPTYPE, typename Device> 
struct resize_memory_op {
  /// @brief Allocate memory for a given pointer. Note this op will free the pointer first.
  ///
  /// Input Parameters
  /// \param dev : the type of computing device
  /// \param size : array size
  ///
  /// Output Parameters
  /// \param arr : allocated array
  void operator()(const Device* dev, FPTYPE*& arr, const size_t size);
};

template <typename FPTYPE, typename Device> 
struct set_memory_op {
  /// @brief memset for multi-device
  ///
  /// Input Parameters
  /// \param dev : the type of computing device
  /// \param var : the specified constant value
  /// \param size : array size
  ///
  /// Output Parameters
  /// \param arr : output array initialized by the input value
  void operator()(const Device* dev, FPTYPE* arr, const int var, const size_t size);
};

template <typename FPTYPE, typename Device_out, typename Device_in> 
struct synchronize_memory_op {
  /// @brief memcpy for multi-device
  ///
  /// Input Parameters
  /// \param dev_out : the type of computing device of arr_out
  /// \param dev_in : the type of computing device of arr_in
  /// \param arr_in : input array
  /// \param size : array size
  ///
  /// Output Parameters
  /// \param arr_out : output array initialized by the input array
  void operator()(
      const Device_out* dev_out, 
      const Device_in* dev_in, 
      FPTYPE* arr_out, 
      const FPTYPE* arr_in, 
      const size_t size);
};

template <typename FPTYPE, typename Device>
struct delete_memory_op {
  /// @brief free memory for multi-device
  ///
  /// Input Parameters
  /// \param dev : the type of computing device
  /// \param arr : the input array
  void operator()(const Device* dev, FPTYPE* arr);
};

#if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
// Partially specialize operator for psi::GpuDevice.
template <typename FPTYPE> 
struct resize_memory_op<FPTYPE, psi::DEVICE_GPU> {
  void operator()(const psi::DEVICE_GPU* dev, FPTYPE*& arr, const size_t size);
};

template <typename FPTYPE> 
struct set_memory_op<FPTYPE, psi::DEVICE_GPU> {
  void operator()(const psi::DEVICE_GPU* dev, FPTYPE* arr, const int var, const size_t size);
};

template <typename FPTYPE> 
struct synchronize_memory_op<FPTYPE, psi::DEVICE_CPU, psi::DEVICE_GPU> {
  void operator()(
      const psi::DEVICE_CPU* dev_out, 
      const psi::DEVICE_GPU* dev_in, 
      FPTYPE* arr_out, 
      const FPTYPE* arr_in, 
      const size_t size);
};
template <typename FPTYPE> 
struct synchronize_memory_op<FPTYPE, psi::DEVICE_GPU, psi::DEVICE_CPU> {
  void operator()(
      const psi::DEVICE_GPU* dev_out, 
      const psi::DEVICE_CPU* dev_in, 
      FPTYPE* arr_out, 
      const FPTYPE* arr_in, 
      const size_t size);
};
template <typename FPTYPE> 
struct synchronize_memory_op<FPTYPE, psi::DEVICE_GPU, psi::DEVICE_GPU> {
  void operator()(
      const psi::DEVICE_GPU* dev_out, 
      const psi::DEVICE_GPU* dev_in, 
      FPTYPE* arr_out, 
      const FPTYPE* arr_in, 
      const size_t size);
};

template <typename FPTYPE>
struct delete_memory_op<FPTYPE, psi::DEVICE_GPU> {
  void operator()(const psi::DEVICE_GPU* dev, FPTYPE* arr);
};
#endif
// __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM 

} // end of namespace memory
} // end of namespace psi

#endif // MODULE_PSI_MEMORY_H_