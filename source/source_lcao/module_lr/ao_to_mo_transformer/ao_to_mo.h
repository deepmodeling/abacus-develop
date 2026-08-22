#pragma once
#include <ATen/core/tensor.h>
#include "source_psi/psi.h"
#include "source_lcao/module_lr/utils/mo_type.h"
#include <vector>
#include <cassert>
#ifdef __MPI
#include "source_base/parallel_2d.h"
#endif
namespace LR
{
    template<typename T>
    void ao_to_mo_forloop_serial(
        const std::vector<container::Tensor>& mat_ao,
        const psi::Psi<T>& coeff,
        const int& nocc,
        const int& nvirt,
        T* const mat_mo,
        const LR_Util::MO_TYPE type = LR_Util::VO);
    template<typename T>
    void ao_to_mo_blas(
        const std::vector<container::Tensor>& mat_ao,
        const psi::Psi<T>& coeff,
        const int& nocc,
        const int& nvirt,
        T* const mat_mo,
        const bool add_on = true,
        const LR_Util::MO_TYPE type = LR_Util::VO);
#ifdef __MPI
    template<typename T>
    void ao_to_mo_pblas(
        const std::vector<container::Tensor>& mat_ao,
        const Parallel_2D& pmat_ao,
        const psi::Psi<T>& coeff,
        const Parallel_2D& pcoeff,
        const int& naos,
        const int& nocc,
        const int& nvirt,
        const Parallel_2D& pmat_mo,
        T* const mat_mo,
        const bool add_on = true,
        const LR_Util::MO_TYPE type = LR_Util::VO);
#endif
}