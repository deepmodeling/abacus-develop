#pragma once
#include <ATen/core/tensor.h>
#include "source_psi/psi.h"
#include <vector>
#ifdef __MPI
#include "source_base/parallel_2d.h"
#endif
namespace LR
{
    // occ
    ///  $\sum_{k\mu\nu}C^*_{\mu i}K_{\mu\nu}C_{\nu k}X_{ak}^*$
    template <typename T>
    void CVCX_occ_forloop_serial(
        const std::vector<container::Tensor>& V_istate,
        const psi::Psi<T, base_device::DEVICE_CPU>& c,
        const T* const X_istate,
        const int& naos,
        const int& nocc,
        const int& nvirt,
        T* const AX_istate);
    template <typename T>
    void CVCX_occ_blas(
        const std::vector<container::Tensor>& V_istate,
        const psi::Psi<T, base_device::DEVICE_CPU>& c,
        const T* const X_istate,
        const int& naos,
        const int& nocc,
        const int& nvirt,
        T* const AX_istate,
        const bool add_on = true,
        const T factor = (T)1.0);
#ifdef __MPI
    template <typename T>
    void CVCX_occ_pblas(
        const std::vector<container::Tensor>& V_istate,
        const Parallel_2D& pmat,
        const psi::Psi<T, base_device::DEVICE_CPU>& c,
        const Parallel_2D& pc,
        const T* const X_istate,
        const Parallel_2D& px,
        const int& naos,
        const int& nocc,
        const int& nvirt,
        T* const AX_istate,
        const bool add_on = true,
        const T factor = (T)1.0);
#endif
    // virt
    ///  $\sum_{b\mu\nu}X^*_{bi}C^*_{\mu b}K_{\mu\nu}C_{\nu a}$
    template <typename T>
    void CVCX_virt_forloop_serial(
        const std::vector<container::Tensor>& V_istate,
        const psi::Psi<T, base_device::DEVICE_CPU>& c,
        const T* const X_istate,
        const int& naos,
        const int& nocc,
        const int& nvirt,
        T* const AX_istate);
    template <typename T>
    void CVCX_virt_blas(
        const std::vector<container::Tensor>& V_istate,
        const psi::Psi<T, base_device::DEVICE_CPU>& c,
        const T* const X_istate,
        const int& naos,
        const int& nocc,
        const int& nvirt,
        T* const AX_istate,
        const bool add_on = true,
        const T factor = (T)1.0);
#ifdef __MPI
    template <typename T>
    void CVCX_virt_pblas(
        const std::vector<container::Tensor>& V_istate,
        const Parallel_2D& pmat,
        const psi::Psi<T, base_device::DEVICE_CPU>& c,
        const Parallel_2D& pc,
        const T* const X_istate,
        const Parallel_2D& px,
        const int& naos,
        const int& nocc,
        const int& nvirt,
        T* const AX_istate,
        const bool add_on = true,
        const T factor = (T)1.0);
#endif
}