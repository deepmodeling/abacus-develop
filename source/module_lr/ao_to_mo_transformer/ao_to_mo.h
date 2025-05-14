#pragma once
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
    template<typename T>
    void  ao_to_mo_forloop_serial(
        const std::vector<container::Tensor>& mat_ao,
        const psi::Psi<T>& coeff,
        const int& nocc,
        const int& nvirt,
        T* const mat_mo,
        const MO_TYPE type = VO);
    template<typename T>
    void ao_to_mo_blas(
        const std::vector<container::Tensor>& mat_ao,
        const psi::Psi<T>& coeff,
        const int& nocc,
        const int& nvirt,
        T* const mat_mo,
        const bool add_on = true,
        const MO_TYPE type = VO);
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
        const MO_TYPE type = VO);
#endif
}