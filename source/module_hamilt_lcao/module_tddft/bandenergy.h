/**
 * @file bandenegy.h
 * @brief compute band energy ekb
 *  This file originally belonged to file LCAO_evolve.cpp
 */
#ifndef BANDENERGY_H
#define BANDENERGY_H

#include "module_base/module_container/ATen/core/tensor.h" // container::Tensor
#include "module_basis/module_ao/parallel_orbitals.h"

#include <complex>

namespace module_tddft
{
#ifdef __MPI
/**
 *  @brief compute band energy ekb <psi_i|H|psi_i>
 *
 * @param[in] pv information of parallel
 * @param[in] nband number of bands
 * @param[in] nlocal number of orbitals
 * @param[in] Htmp Hamiltonian
 * @param[in] psi_k psi of this step
 * @param[out] ekb band energy
 */
void compute_ekb(const Parallel_Orbitals* pv,
                 const int nband,
                 const int nlocal,
                 const std::complex<double>* Htmp,
                 const std::complex<double>* psi_k,
                 double* ekb);

void compute_ekb_tensor(const Parallel_Orbitals* pv,
                        const int nband,
                        const int nlocal,
                        const container::Tensor& Htmp,
                        const container::Tensor& psi_k,
                        container::Tensor& ekb);

void compute_ekb_tensor_lapack(const Parallel_Orbitals* pv,
                               const int nband,
                               const int nlocal,
                               const container::Tensor& Htmp,
                               const container::Tensor& psi_k,
                               container::Tensor& ekb);
#endif
} // namespace module_tddft
#endif
