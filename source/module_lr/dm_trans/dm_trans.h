#pragma once
// use tensor or basematrix in the future
#include <ATen/core/tensor.h>
#include "module_psi/psi.h"
#include <vector>
#ifdef __MPI
#include "module_base/parallel_2d.h"
#endif
namespace LR
{

#ifndef MO_TYPE_H
#define MO_TYPE_H
    enum MO_TYPE { OO, VO, VV };
#endif

#ifdef __MPI
/// @brief calculate the 2d-block transition density matrix in AO basis using p?gemm
/// \f[ \tilde{\rho}_{\mu_j\mu_b}=\sum_{jb}c_{j,\mu_j}X_{jb}c^*_{b,\mu_b} \f]
    template<typename T>
    std::vector<container::Tensor> cal_dm_trans_pblas(
        const T* const X_istate,
        const Parallel_2D& px,
        const psi::Psi<T>& c,
        const Parallel_2D& pc,
        const int naos,
        const int nocc,
        const int nvirt,
        const Parallel_2D& pmat,
        const T factor = (T)1.0,
        const MO_TYPE type = MO_TYPE::VO);
#endif

    /// @brief calculate the 2d-block transition density matrix in AO basis using ?gemm
    template<typename T>
    std::vector<container::Tensor> cal_dm_trans_blas(
        const T* const X_istate,
        const psi::Psi<T>& c,
        const int& nocc, const int& nvirt,
        const T factor = (T)1.0,
        const MO_TYPE type = MO_TYPE::VO);

    // for test
    /// @brief calculate the 2d-block transition density matrix in AO basis using for loop (for test)
    template<typename T>
    std::vector<container::Tensor> cal_dm_trans_forloop_serial(
        const T* const X_istate,
        const psi::Psi<T>& c,
        const int& nocc, const int& nvirt,
        const T factor = (T)1.0,
        const MO_TYPE type = MO_TYPE::VO);
}