#ifndef MODULE_DEVICE_MEMORY_H_
#define MODULE_DEVICE_MEMORY_H_

#include "device.h"

#include <complex>
#include <stddef.h>
#include <vector>

namespace device
{

namespace memory
{

template <typename FPTYPE, typename Device>
struct resize_memory_op
{
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
struct set_memory_op
{
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
struct synchronize_memory_op
{
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
    void operator()(const Device_out* dev_out,
                    const Device_in* dev_in,
                    FPTYPE* arr_out,
                    const FPTYPE* arr_in,
                    const size_t size);
};

template <typename FPTYPE_out, typename FPTYPE_in, typename Device_out, typename Device_in>
struct cast_memory_op
{
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
    void operator()(const Device_out* dev_out,
                    const Device_in* dev_in,
                    FPTYPE_out* arr_out,
                    const FPTYPE_in* arr_in,
                    const size_t size);
};

template <typename FPTYPE, typename Device>
struct delete_memory_op
{
    /// @brief free memory for multi-device
    ///
    /// Input Parameters
    /// \param dev : the type of computing device
    /// \param arr : the input array
    void operator()(const Device* dev, FPTYPE* arr);
};
} // end of namespace memory
} // end of namespace device

using resmem_sh_op = device::memory::resize_memory_op<float, device::CPU>;
using resmem_dh_op = device::memory::resize_memory_op<double, device::CPU>;
using resmem_ch_op = device::memory::resize_memory_op<std::complex<float>, device::CPU>;
using resmem_zh_op = device::memory::resize_memory_op<std::complex<double>, device::CPU>;

using resmem_sd_op = device::memory::resize_memory_op<float, device::GPU>;
using resmem_dd_op = device::memory::resize_memory_op<double, device::GPU>;
using resmem_cd_op = device::memory::resize_memory_op<std::complex<float>, device::GPU>;
using resmem_zd_op = device::memory::resize_memory_op<std::complex<double>, device::GPU>;

using setmem_sh_op = device::memory::set_memory_op<float, device::CPU>;
using setmem_dh_op = device::memory::set_memory_op<double, device::CPU>;
using setmem_ch_op = device::memory::set_memory_op<std::complex<float>, device::CPU>;
using setmem_zh_op = device::memory::set_memory_op<std::complex<double>, device::CPU>;

using setmem_sd_op = device::memory::set_memory_op<float, device::GPU>;
using setmem_dd_op = device::memory::set_memory_op<double, device::GPU>;
using setmem_cd_op = device::memory::set_memory_op<std::complex<float>, device::GPU>;
using setmem_zd_op = device::memory::set_memory_op<std::complex<double>, device::GPU>;

using delmem_sh_op = device::memory::delete_memory_op<float, device::CPU>;
using delmem_dh_op = device::memory::delete_memory_op<double, device::CPU>;
using delmem_ch_op = device::memory::delete_memory_op<std::complex<float>, device::CPU>;
using delmem_zh_op = device::memory::delete_memory_op<std::complex<double>, device::CPU>;

using delmem_sd_op = device::memory::delete_memory_op<float, device::GPU>;
using delmem_dd_op = device::memory::delete_memory_op<double, device::GPU>;
using delmem_cd_op = device::memory::delete_memory_op<std::complex<float>, device::GPU>;
using delmem_zd_op = device::memory::delete_memory_op<std::complex<double>, device::GPU>;

using syncmem_s2s_h2h_op = device::memory::synchronize_memory_op<float, device::CPU, device::CPU>;
using syncmem_s2s_h2d_op = device::memory::synchronize_memory_op<float, device::GPU, device::CPU>;
using syncmem_s2s_d2h_op = device::memory::synchronize_memory_op<float, device::CPU, device::GPU>;
using syncmem_d2d_h2h_op = device::memory::synchronize_memory_op<double, device::CPU, device::CPU>;
using syncmem_d2d_h2d_op = device::memory::synchronize_memory_op<double, device::GPU, device::CPU>;
using syncmem_d2d_d2h_op = device::memory::synchronize_memory_op<double, device::CPU, device::GPU>;

using syncmem_c2c_h2h_op
    = device::memory::synchronize_memory_op<std::complex<float>, device::CPU, device::CPU>;
using syncmem_c2c_h2d_op
    = device::memory::synchronize_memory_op<std::complex<float>, device::GPU, device::CPU>;
using syncmem_c2c_d2h_op
    = device::memory::synchronize_memory_op<std::complex<float>, device::CPU, device::GPU>;
using syncmem_z2z_h2h_op
    = device::memory::synchronize_memory_op<std::complex<double>, device::CPU, device::CPU>;
using syncmem_z2z_h2d_op
    = device::memory::synchronize_memory_op<std::complex<double>, device::GPU, device::CPU>;
using syncmem_z2z_d2h_op
    = device::memory::synchronize_memory_op<std::complex<double>, device::CPU, device::GPU>;

using castmem_s2d_h2h_op = device::memory::cast_memory_op<double, float, device::CPU, device::CPU>;
using castmem_s2d_h2d_op = device::memory::cast_memory_op<double, float, device::GPU, device::CPU>;
using castmem_s2d_d2h_op = device::memory::cast_memory_op<double, float, device::CPU, device::GPU>;
using castmem_d2s_h2h_op = device::memory::cast_memory_op<float, double, device::CPU, device::CPU>;
using castmem_d2s_h2d_op = device::memory::cast_memory_op<float, double, device::GPU, device::CPU>;
using castmem_d2s_d2h_op = device::memory::cast_memory_op<float, double, device::CPU, device::GPU>;

using castmem_c2z_h2h_op
    = device::memory::cast_memory_op<std::complex<double>, std::complex<float>, device::CPU, device::CPU>;
using castmem_c2z_h2d_op
    = device::memory::cast_memory_op<std::complex<double>, std::complex<float>, device::GPU, device::CPU>;
using castmem_c2z_d2h_op
    = device::memory::cast_memory_op<std::complex<double>, std::complex<float>, device::CPU, device::GPU>;
using castmem_z2c_h2h_op
    = device::memory::cast_memory_op<std::complex<float>, std::complex<double>, device::CPU, device::CPU>;
using castmem_z2c_h2d_op
    = device::memory::cast_memory_op<std::complex<float>, std::complex<double>, device::GPU, device::CPU>;
using castmem_z2c_d2h_op
    = device::memory::cast_memory_op<std::complex<float>, std::complex<double>, device::CPU, device::GPU>;

static device::CPU* cpu_ctx = {};
static device::GPU* gpu_ctx = {};

#endif // MODULE_DEVICE_MEMORY_H_