#ifndef NORM_PSI_H
#define NORM_PSI_H

#include "module_basis/module_ao/parallel_orbitals.h"

namespace module_tddft
{
#ifdef __MPI
void norm_psi(const Parallel_Orbitals* pv,
              const int nband,
              const int nlocal,
              const std::complex<double>* Stmp,
              std::complex<double>* psi_k,
              const int print_matrix);

#endif
} // namespace module_tddft

#endif