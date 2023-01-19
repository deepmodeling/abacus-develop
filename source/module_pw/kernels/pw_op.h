#ifndef MODULE_PW_MULTI_DEVICE_H
#define MODULE_PW_MULTI_DEVICE_H

#include "module_psi/psi.h"
#include <complex>

namespace ModulePW {
template <typename FPTYPE, typename Device>
struct set_3d_fft_box_op {
    /// @brief Set the 3D fft box for fft transfrom between the recip and real space.
    /// To map the 1D psi(1D continuous array) to 3D box psi(fft box)
    ///
    /// Input Parameters
    /// @param dev - which device this function runs on
    /// @param npwk - number of planwaves
    /// @param box_index - the mapping function of 1D to 3D
    /// @param in - input psi within a 1D array(in recip space)
    ///
    /// Output Parameters
    /// @param out - output psi within the 3D box(in recip space)
    void operator() (
        const Device* dev,
        const int npwk,
        const int* box_index,
        const std::complex<FPTYPE>* in,
        std::complex<FPTYPE>* out);
};

template <typename FPTYPE, typename Device>
struct set_recip_to_real_output_op {
    /// @brief Calculate the outputs after the FFT translation of recip_to_real
    ///
    /// Input Parameters
    /// @param dev - which device this function runs on
    /// @param nrxx - size of array
    /// @param add - flag to control whether to add the input itself
    /// @param in - input psi within a 1D array(in real space)
    ///
    /// Output Parameters
    /// @param out - output psi within the 3D box(in real space)
    void operator() (
        const Device* dev,
        const int nrxx,
        const bool add,
        const FPTYPE factor,
        const std::complex<FPTYPE>* in,
        std::complex<FPTYPE>* out);
};

template <typename FPTYPE, typename Device>
struct set_real_to_recip_output_op {
    /// @brief Calculate the outputs after the FFT translation of real_to_recip
    ///
    /// Input Parameters
    /// @param dev - which device this function runs on
    /// @param nxyz - size of array
    /// @param add - flag to control whether to add the input itself
    /// @param factor - input constant value
    /// @param box_index - input box parameters
    /// @param in - input psi within a 1D array(in recip space)
    ///
    /// Output Parameters
    /// @param out - output psi within the 3D box(in recip space)
    void operator() (
        const Device* dev,
        const int npw_k,
        const int nxyz,
        const bool add,
        const FPTYPE factor,
        const int* box_index,
        const std::complex<FPTYPE>* in,
        std::complex<FPTYPE>* out);
};

#if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
// Partially specialize functor for psi::GpuDevice.
template <typename FPTYPE>
struct set_3d_fft_box_op<FPTYPE, psi::DEVICE_GPU> {
    void operator() (
        const psi::DEVICE_GPU* dev,
        const int npwk,
        const int* box_index,
        const std::complex<FPTYPE>* in,
        std::complex<FPTYPE>* out);
};

template <typename FPTYPE>
struct set_recip_to_real_output_op<FPTYPE, psi::DEVICE_GPU> {
    void operator() (
        const psi::DEVICE_GPU* dev,
        const int nrxx,
        const bool add,
        const FPTYPE factor,
        const std::complex<FPTYPE>* in,
        std::complex<FPTYPE>* out);
};

template <typename FPTYPE>
struct set_real_to_recip_output_op<FPTYPE, psi::DEVICE_GPU> {
    void operator() (
        const psi::DEVICE_GPU* dev,
        const int npw_k,
        const int nxyz,
        const bool add,
        const FPTYPE factor,
        const int* box_index,
        const std::complex<FPTYPE>* in,
        std::complex<FPTYPE>* out);
};

#endif // __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
} // namespace ModulePW
#endif //MODULE_PW_MULTI_DEVICE_H