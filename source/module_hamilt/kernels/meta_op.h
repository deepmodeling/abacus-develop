#ifndef MODULE_HAMILT_META_H
#define MODULE_HAMILT_META_H

#include "module_psi/psi.h"
#include <complex>

namespace hamilt {
    template <typename FPTYPE, typename Device>
    struct meta_pw_op {
        /// @brief Compute the metaGGA potential of hPsi in recip space,
        /// out = hpsi * (wfcpw->gcar+wfcpw->kvec_c) * this->tpiba;
        ///
        /// Input Parameters
        /// \param dev : the type of computing device
        /// \param ik : current k point
        /// \param pol : loop of [0, 1, 2]
        /// \param npw : number of planewaves of current k point
        /// \param npwx : number of planewaves of all k points
        /// \param tpiba : GlobalC::ucell.tpiba
        /// \param gcar : wfcpw->gcar
        /// \param kvec_c : wfcpw->kvec_c
        /// \param in : input hpsi
        /// \param bool : conditional flag
        ///
        /// Output Parameters
        /// \param out : output array with size of nspin * npwx
        void operator() (
                const Device* dev,
                const int& ik,
                const int& pol,
                const int& npw,
                const int& npwx,
                const FPTYPE& tpiba,
                const FPTYPE* gcar,
                const FPTYPE* kvec_c,
                const std::complex<FPTYPE>* in,
                std::complex<FPTYPE>* out,
                const bool add = false);
    };

#if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
// Partially specialize functor for psi::GpuDevice.
    template <typename FPTYPE>
    struct meta_pw_op<FPTYPE, psi::DEVICE_GPU> {
        void operator() (
                const psi::DEVICE_GPU* dev,
                const int& ik,
                const int& pol,
                const int& npw,
                const int& npwx,
                const FPTYPE& tpiba,
                const FPTYPE* gcar,
                const FPTYPE* kvec_c,
                const std::complex<FPTYPE>* in,
                std::complex<FPTYPE>* out,
                const bool add = false);
    };
#endif // __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
} // namespace hamilt
#endif //MODULE_HAMILT_META_H