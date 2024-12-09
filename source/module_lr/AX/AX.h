#pragma once
#include <ATen/core/tensor.h>
#include "module_psi/psi.h"
#include <vector>
#ifdef __MPI
#include "module_base/parallel_2d.h"
#endif
namespace LR
{
    // double
    void  cal_AX_forloop_serial(
        const std::vector<container::Tensor>& V_istate,
        const psi::Psi<double>& c,
        const int& nocc,
        const int& nvirt,
        double* AX_istate);
    void cal_AX_blas(
        const std::vector<container::Tensor>& V_istate,
        const psi::Psi<double>& c,
        const int& nocc,
        const int& nvirt,
        double* AX_istate,
        const bool add_on = true);
#ifdef __MPI
    void cal_AX_pblas(
        const std::vector<container::Tensor>& V_istate,
        const Parallel_2D& pmat,
        const psi::Psi<double>& c,
        const Parallel_2D& pc,
        const int& naos,
        const int& nocc,
        const int& nvirt,
        const Parallel_2D& pX,
        double* AX_istate,
        const bool add_on=true);
#endif
    // complex
    void cal_AX_forloop_serial(
        const std::vector<container::Tensor>& V_istate,
        const psi::Psi<std::complex<double>>& c,
        const int& nocc,
        const int& nvirt,
        std::complex<double>* AX_istate);
    void cal_AX_blas(
        const std::vector<container::Tensor>& V_istate,
        const psi::Psi<std::complex<double>>& c,
        const int& nocc,
        const int& nvirt,
        std::complex<double>* AX_istate,
        const bool add_on = true);

#ifdef __MPI
    void  cal_AX_pblas(
        const std::vector<container::Tensor>& V_istate,
        const Parallel_2D& pmat,
        const psi::Psi<std::complex<double>>& c,
        const Parallel_2D& pc,
        const int& naos,
        const int& nocc,
        const int& nvirt,
        const Parallel_2D& pX,
        std::complex<double>* AX_istate,
        const bool add_on = true);
#endif
}