// TODO: This is a temperary location for these functions.
// And will be moved to a global module(module base) later.
#ifndef MODULE_PSI_MEMORY_H_
#define MODULE_PSI_MEMORY_H_

#include "types.h"

#include <vector>
#include <complex>
#include <stddef.h>

namespace psi {
namespace memory {

template <typename FPTYPE, typename Device> 
struct resize_memory_op {
  /// @brief Allocate memory for a given pointer. Note this op will free the pointer first.
  ///
  /// Input Parameters
  /// \param dev : the type of computing device
  /// \param size : array size
  /// \param record_string : label for memory record
  ///
  /// Output Parameters
  /// \param arr : allocated array
  void operator()(const Device* dev, FPTYPE*& arr, const size_t size, const char* record_in = nullptr);
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

template <typename FPTYPE_out, typename FPTYPE_in, typename Device_out, typename Device_in>
struct cast_memory_op {
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
            FPTYPE_out* arr_out,
            const FPTYPE_in* arr_in,
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
  void operator()(const psi::DEVICE_GPU* dev, FPTYPE*& arr, const size_t size, const char* record_in = nullptr);
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

using resmem_sh_op = psi::memory::resize_memory_op<float, psi::DEVICE_CPU>;
using resmem_dh_op = psi::memory::resize_memory_op<double, psi::DEVICE_CPU>;
using resmem_ch_op = psi::memory::resize_memory_op<std::complex<float>, psi::DEVICE_CPU>;
using resmem_zh_op = psi::memory::resize_memory_op<std::complex<double>, psi::DEVICE_CPU>;

using resmem_sd_op = psi::memory::resize_memory_op<float, psi::DEVICE_GPU>;
using resmem_dd_op = psi::memory::resize_memory_op<double, psi::DEVICE_GPU>;
using resmem_cd_op = psi::memory::resize_memory_op<std::complex<float>, psi::DEVICE_GPU>;
using resmem_zd_op = psi::memory::resize_memory_op<std::complex<double>, psi::DEVICE_GPU>;

using setmem_sh_op = psi::memory::set_memory_op<float, psi::DEVICE_CPU>;
using setmem_dh_op = psi::memory::set_memory_op<double, psi::DEVICE_CPU>;
using setmem_ch_op = psi::memory::set_memory_op<std::complex<float>, psi::DEVICE_CPU>;
using setmem_zh_op = psi::memory::set_memory_op<std::complex<double>, psi::DEVICE_CPU>;

using setmem_sd_op = psi::memory::set_memory_op<float, psi::DEVICE_GPU>;
using setmem_dd_op = psi::memory::set_memory_op<double, psi::DEVICE_GPU>;
using setmem_cd_op = psi::memory::set_memory_op<std::complex<float>, psi::DEVICE_GPU>;
using setmem_zd_op = psi::memory::set_memory_op<std::complex<double>, psi::DEVICE_GPU>;

using delmem_sh_op = psi::memory::delete_memory_op<float, psi::DEVICE_CPU>;
using delmem_dh_op = psi::memory::delete_memory_op<double, psi::DEVICE_CPU>;
using delmem_ch_op = psi::memory::delete_memory_op<std::complex<float>, psi::DEVICE_CPU>;
using delmem_zh_op = psi::memory::delete_memory_op<std::complex<double>, psi::DEVICE_CPU>;

using delmem_sd_op = psi::memory::delete_memory_op<float, psi::DEVICE_GPU>;
using delmem_dd_op = psi::memory::delete_memory_op<double, psi::DEVICE_GPU>;
using delmem_cd_op = psi::memory::delete_memory_op<std::complex<float>, psi::DEVICE_GPU>;
using delmem_zd_op = psi::memory::delete_memory_op<std::complex<double>, psi::DEVICE_GPU>;

using syncmem_s2s_h2h_op = psi::memory::synchronize_memory_op<float, psi::DEVICE_CPU, psi::DEVICE_CPU>;
using syncmem_s2s_h2d_op = psi::memory::synchronize_memory_op<float, psi::DEVICE_GPU, psi::DEVICE_CPU>;
using syncmem_s2s_d2h_op = psi::memory::synchronize_memory_op<float, psi::DEVICE_CPU, psi::DEVICE_GPU>;
using syncmem_d2d_h2h_op = psi::memory::synchronize_memory_op<double, psi::DEVICE_CPU, psi::DEVICE_CPU>;
using syncmem_d2d_h2d_op = psi::memory::synchronize_memory_op<double, psi::DEVICE_GPU, psi::DEVICE_CPU>;
using syncmem_d2d_d2h_op = psi::memory::synchronize_memory_op<double, psi::DEVICE_CPU, psi::DEVICE_GPU>;

using syncmem_c2c_h2h_op = psi::memory::synchronize_memory_op<std::complex<float>, psi::DEVICE_CPU, psi::DEVICE_CPU>;
using syncmem_c2c_h2d_op = psi::memory::synchronize_memory_op<std::complex<float>, psi::DEVICE_GPU, psi::DEVICE_CPU>;
using syncmem_c2c_d2h_op = psi::memory::synchronize_memory_op<std::complex<float>, psi::DEVICE_CPU, psi::DEVICE_GPU>;
using syncmem_z2z_h2h_op = psi::memory::synchronize_memory_op<std::complex<double>, psi::DEVICE_CPU, psi::DEVICE_CPU>;
using syncmem_z2z_h2d_op = psi::memory::synchronize_memory_op<std::complex<double>, psi::DEVICE_GPU, psi::DEVICE_CPU>;
using syncmem_z2z_d2h_op = psi::memory::synchronize_memory_op<std::complex<double>, psi::DEVICE_CPU, psi::DEVICE_GPU>;

using castmem_s2d_h2h_op = psi::memory::cast_memory_op<double, float, psi::DEVICE_CPU, psi::DEVICE_CPU>;
using castmem_s2d_h2d_op = psi::memory::cast_memory_op<double, float, psi::DEVICE_GPU, psi::DEVICE_CPU>;
using castmem_s2d_d2h_op = psi::memory::cast_memory_op<double, float, psi::DEVICE_CPU, psi::DEVICE_GPU>;
using castmem_d2s_h2h_op = psi::memory::cast_memory_op<float, double, psi::DEVICE_CPU, psi::DEVICE_CPU>;
using castmem_d2s_h2d_op = psi::memory::cast_memory_op<float, double, psi::DEVICE_GPU, psi::DEVICE_CPU>;
using castmem_d2s_d2h_op = psi::memory::cast_memory_op<float, double, psi::DEVICE_CPU, psi::DEVICE_GPU>;

using castmem_c2z_h2h_op = psi::memory::cast_memory_op<std::complex<double>, std::complex<float>, psi::DEVICE_CPU, psi::DEVICE_CPU>;
using castmem_c2z_h2d_op = psi::memory::cast_memory_op<std::complex<double>, std::complex<float>, psi::DEVICE_GPU, psi::DEVICE_CPU>;
using castmem_c2z_d2h_op = psi::memory::cast_memory_op<std::complex<double>, std::complex<float>, psi::DEVICE_CPU, psi::DEVICE_GPU>;
using castmem_z2c_h2h_op = psi::memory::cast_memory_op<std::complex<float>, std::complex<double>, psi::DEVICE_CPU, psi::DEVICE_CPU>;
using castmem_z2c_h2d_op = psi::memory::cast_memory_op<std::complex<float>, std::complex<double>, psi::DEVICE_GPU, psi::DEVICE_CPU>;
using castmem_z2c_d2h_op = psi::memory::cast_memory_op<std::complex<float>, std::complex<double>, psi::DEVICE_CPU, psi::DEVICE_GPU>;

static psi::DEVICE_CPU * cpu_ctx = {};
static psi::DEVICE_GPU * gpu_ctx = {};

#endif // MODULE_PSI_MEMORY_H_