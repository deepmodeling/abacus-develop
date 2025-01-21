/**
 * @file middle_hamilt.h
 * @brief  compute H(t+dt/2)
 *  This file originally belonged to file LCAO_evolve.cpp
 */
#ifndef MIDDLE_HAMILT_H
#define MIDDLE_HAMILT_H

#include "module_base/module_container/ATen/core/tensor.h" // ct::Tensor
#include "module_basis/module_ao/parallel_orbitals.h"

#include <complex>

namespace module_tddft
{
#ifdef __MPI
/**
 *  @brief compute H(t+dt/2)
 *
 * @param[in] pv information of parallel
 * @param[in] nband number of bands
 * @param[in] nlocal number of orbitals
 * @param[in] Htmp H(t+dt)
 * @param[in] H_laststep H(t)
 * @param[in] print_matirx print internal matrix or not
 * @param[out] Htmp H(t+dt/2)
 */
void half_Hmatrix(const Parallel_Orbitals* pv,
                  const int nband,
                  const int nlocal,
                  std::complex<double>* Htmp,
                  std::complex<double>* Stmp,
                  const std::complex<double>* H_laststep,
                  const std::complex<double>* S_laststep,
                  std::ofstream& ofs_running,
                  const int print_matrix);

void half_Hmatrix_tensor(const Parallel_Orbitals* pv,
                         const int nband,
                         const int nlocal,
                         ct::Tensor& Htmp,
                         ct::Tensor& Stmp,
                         const ct::Tensor& H_laststep,
                         const ct::Tensor& S_laststep,
                         std::ofstream& ofs_running,
                         const int print_matrix);

template <typename Device>
void half_Hmatrix_tensor_lapack(const Parallel_Orbitals* pv,
                                const int nband,
                                const int nlocal,
                                ct::Tensor& Htmp,
                                ct::Tensor& Stmp,
                                const ct::Tensor& H_laststep,
                                const ct::Tensor& S_laststep,
                                std::ofstream& ofs_running,
                                const int print_matrix);

#endif // __MPI
} // namespace module_tddft

#endif
