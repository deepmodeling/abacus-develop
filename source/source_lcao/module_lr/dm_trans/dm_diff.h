#ifndef ABACUS_SOURCE_LCAO_MODULE_LR_DM_TRANS_DM_DIFF_H
#define ABACUS_SOURCE_LCAO_MODULE_LR_DM_TRANS_DM_DIFF_H
#include <ATen/core/tensor.h>
#include "source_psi/psi.h"
#include <vector>
#include "source_lcao/module_lr/utils/lr_util.h"
namespace LR
{
    /// Orbital-amplitude contractions C_v X and C_o X^T, shared with the EXX operator.
    void CvX(
        const double* C,
        const double* X,
        const int& naos,
        const int& nocc,
        const int& nvirt,
        double* CvX);
    void CvX(
        const std::complex<double>* C,
        const std::complex<double>* X,
        const int& naos,
        const int& nocc,
        const int& nvirt,
        std::complex<double>* CvX);
    void CoXT(
        const double* C,
        const double* X,
        const int& naos,
        const int& nocc,
        const int& nvirt,
        double* CoXT);
    void CoXT(
        const std::complex<double>* C,
        const std::complex<double>* X,
        const int& naos,
        const int& nocc,
        const int& nvirt,
        std::complex<double>* CoXT);
#ifdef __MPI
    void CvX(
        const double* C,
        const BlockCyclicDistributionInformation& pc,
        const double* X,
        const BlockCyclicDistributionInformation& px,
        const int& naos,
        const int& nocc,
        const int& nvirt,
        double* CvX,
        const BlockCyclicDistributionInformation& pcx);
    void CvX(
        const std::complex<double>* C,
        const BlockCyclicDistributionInformation& pc,
        const std::complex<double>* X,
        const BlockCyclicDistributionInformation& px,
        const int& naos,
        const int& nocc,
        const int& nvirt,
        std::complex<double>* CvX,
        const BlockCyclicDistributionInformation& pcx);
    void CoXT(
        const double* C,
        const BlockCyclicDistributionInformation& pc,
        const double* X,
        const BlockCyclicDistributionInformation& px,
        const int& naos,
        const int& nocc,
        const int& nvirt,
        double* CoXT,
        const BlockCyclicDistributionInformation& pcxt);
    void CoXT(
        const std::complex<double>* C,
        const BlockCyclicDistributionInformation& pc,
        const std::complex<double>* X,
        const BlockCyclicDistributionInformation& px,
        const int& naos,
        const int& nocc,
        const int& nvirt,
        std::complex<double>* CoXT,
        const BlockCyclicDistributionInformation& pcxt);
#endif
    // use templates in the future.
#ifdef __MPI
/// @brief calculate the 2d-block difference density matrix in AO basis using p?gemm
/// \f[ T=(C_v*X) * (C_v * X)^\dagger + (C_o*X^T) * (C_o*X^T)^\dagger \f]
    template<typename T>
    std::vector<container::Tensor> cal_dm_diff_pblas(
        const T* const X_istate,
        const BlockCyclicDistributionInformation& px,
        const psi::Psi<T, base_device::DEVICE_CPU>& c,
        const BlockCyclicDistributionInformation& pc,
        const int& naos,
        const int& nocc,
        const int& nvirt,
        const BlockCyclicDistributionInformation& pmat,
        const bool renorm_k = true,
        const int nspin = 1);
#endif

    /// @brief calculate the 2d-block transition density matrix in AO basis using ?gemm
    template<typename T>
    std::vector<container::Tensor> cal_dm_diff_blas(
        const T* const X_istate,
        const psi::Psi<T, base_device::DEVICE_CPU>& c,
        const int& naos,
        const int& nocc,
        const int& nvirt,
        const bool renorm_k = true,
        const int nspin = 1);

    // for test
    /// @brief calculate the 2d-block transition density matrix in AO basis using for loop (for test)
    template<typename T>
    std::vector<container::Tensor> cal_dm_diff_forloop(
        const T* const X_istate,
        const psi::Psi<T, base_device::DEVICE_CPU>& c,
        const int& naos,
        const int& nocc,
        const int& nvirt,
        const bool renorm_k = true,
        const int nspin = 1);
}

#endif // ABACUS_SOURCE_LCAO_MODULE_LR_DM_TRANS_DM_DIFF_H
