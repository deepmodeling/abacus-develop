/**
 * @file upsi.h
 * @brief apply U_operator to the wave function of the previous step for new wave function
 *  This file originally belonged to file LCAO_evolve.cpp
 */

#ifndef UPSI_H
#define UPSI_H

#include "module_basis/module_ao/parallel_orbitals.h"

namespace module_tddft
{
#ifdef __MPI
/**
 *  @brief apply U_operator to the wave function of the previous step for new wave function
 *
 * @param[in] pv information of parallel
 * @param[in] nband number of bands
 * @param[in] nlocal number of orbitals
 * @param[in] U_operator operator of propagator
 * @param[in] psi_k_laststep psi of last step
 * @param[in] print_matirx print internal matrix or not
 * @param[out] psi_k psi of this step
 */
void upsi(const Parallel_Orbitals* pv,
          const int nband,
          const int nlocal,
          const std::complex<double>* U_operator,
          const std::complex<double>* psi_k_laststep,
          std::complex<double>* psi_k,
          const int print_matrix);

#endif
} // namespace module_tddft

#endif