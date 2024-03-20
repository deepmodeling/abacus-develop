/**
 * @file norm_psi.h
 * @brief normalize the wave function
 *  This file originally belonged to file LCAO_evolve.cpp
 */
#ifndef NORM_PSI_H
#define NORM_PSI_H

#include "module_basis/module_ao/parallel_orbitals.h"

namespace module_tddft
{
#ifdef __MPI
/**
 *  @brief normalize the wave function
 *
 * @param[in] pv information of parallel
 * @param[in] nband number of bands
 * @param[in] nlocal number of orbitals
 * @param[in] Stmp overlap matrix
 * @param[in] psi_k psi of this step
 * @param[in] print_matirx print internal matrix or not
 * @param[out] psi_k psi of this step after normalization
 */
void norm_psi(const Parallel_Orbitals* pv,
              const int nband,
              const int nlocal,
              const std::complex<double>* Stmp,
              std::complex<double>* psi_k,
              const int print_matrix);

#endif
} // namespace module_tddft

#endif