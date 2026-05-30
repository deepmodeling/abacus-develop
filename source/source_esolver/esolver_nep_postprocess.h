#ifndef ESOLVER_NEP_POSTPROCESS_H
#define ESOLVER_NEP_POSTPROCESS_H

#include "source_base/matrix.h"

namespace ModuleESolver
{

void postprocess_nep_cpu(const int nat,
                         const double* atomic_energy,
                         const double* raw_force,
                         const double* raw_virial,
                         const double fact_e,
                         const double fact_f,
                         const double fact_v,
                         double& potential,
                         ModuleBase::matrix& force,
                         ModuleBase::matrix& virial);

#ifdef __CUDA
void postprocess_nep_cuda(const int nat,
                          const double* atomic_energy,
                          const double* raw_force,
                          const double* raw_virial,
                          const double fact_e,
                          const double fact_f,
                          const double fact_v,
                          double& potential,
                          ModuleBase::matrix& force,
                          ModuleBase::matrix& virial);
#endif

} // namespace ModuleESolver

#endif // ESOLVER_NEP_POSTPROCESS_H
