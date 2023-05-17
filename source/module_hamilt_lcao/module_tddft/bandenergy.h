/**
 * @file bandenegy.h
 * @brief compute band energy ekb
 *  This file originally belonged to file LCAO_evolve.cpp
 */
#ifndef BANDENERGY_H
#define BANDENERGY_H

#include "module_basis/module_ao/parallel_orbitals.h"

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
#endif
} // namespace module_tddft
#endif