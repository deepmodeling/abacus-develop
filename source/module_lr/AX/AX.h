#pragma once
#include <ATen/core/tensor.h>
#include "module_psi/psi.h"
#include <vector>
#ifdef __MPI
#include "module_basis/module_ao/parallel_2d.h"
#endif
namespace LR
{
    // double
    void  cal_AX_forloop_serial(
        const std::vector<container::Tensor>& V_istate,
        const psi::Psi<double>& c,
        const int& nocc,
        const int& nvirt,
        psi::Psi<double>& AX_istate);
    void cal_AX_blas(
        const std::vector<container::Tensor>& V_istate,
        const psi::Psi<double>& c,
        const int& nocc,
        const int& nvirt,
        psi::Psi<double>& AX_istate,
        const bool add_on = true);
#ifdef __MPI
    void cal_AX_pblas(
        const std::vector<container::Tensor>& V_istate,
        const Parallel_2D& pmat,
        const psi::Psi<double>& c,
        const Parallel_2D& pc,
        int naos,
        int nocc,
        int nvirt,
        Parallel_2D& pX,
        psi::Psi<double>& AX_istate,
        const bool add_on=true);
#endif
    // complex
    void cal_AX_forloop_serial(
        const std::vector<container::Tensor>& V_istate,
        const psi::Psi<std::complex<double>>& c,
        const int& nocc,
        const int& nvirt,
        psi::Psi<std::complex<double>>& AX_istate);
    void cal_AX_blas(
        const std::vector<container::Tensor>& V_istate,
        const psi::Psi<std::complex<double>>& c,
        const int& nocc,
        const int& nvirt,
        psi::Psi<std::complex<double>>& AX_istate,
        const bool add_on = true);

#ifdef __MPI
    void  cal_AX_pblas(
        const std::vector<container::Tensor>& V_istate,
        const Parallel_2D& pmat,
        const psi::Psi<std::complex<double>>& c,
        const Parallel_2D& pc,
        int naos,
        int nocc,
        int nvirt,
        Parallel_2D& pX,
        psi::Psi<std::complex<double>>& AX_istate,
        const bool add_on = true);
#endif
}