#ifndef BANDENERGY_H
#define BANDENERGY_H

#include "module_basis/module_ao/parallel_orbitals.h"

#ifdef __MPI
void compute_ekb(const Parallel_Orbitals* pv,
                 const int nband,
                 const int nlocal,
                 const std::complex<double>* Htmp,
                 const std::complex<double>* psi_k,
                 double* ekb);
#endif

#endif