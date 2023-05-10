#ifndef MIDDLE_HAMILT_H
#define MIDDLE_HAMILT_H

#include "module_basis/module_ao/parallel_orbitals.h"

#ifdef __MPI
void half_Hmatrix(const Parallel_Orbitals* pv,
                  const int nband,
                  const int nlocal,
                  std::complex<double>* Htmp,
                  const std::complex<double>* H_laststep,
                  const int print_matrix);
#endif

#endif