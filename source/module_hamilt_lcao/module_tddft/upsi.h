#ifndef UPSI_H
#define UPSI_H

#include "module_basis/module_ao/parallel_orbitals.h"

namespace module_tddft
{
#ifdef __MPI
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