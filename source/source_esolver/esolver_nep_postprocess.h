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
struct NepCudaPostprocessWorkspace
{
    int capacity = 0;
    double* energy = nullptr;
    double* raw_force = nullptr;
    double* raw_virial = nullptr;
    double* potential = nullptr;
    double* force = nullptr;
    double* virial = nullptr;
};

void init_nep_cuda_postprocess_workspace(NepCudaPostprocessWorkspace& workspace, const int nat);

void release_nep_cuda_postprocess_workspace(NepCudaPostprocessWorkspace& workspace);

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

void postprocess_nep_cuda(const int nat,
                          const double* atomic_energy,
                          const double* raw_force,
                          const double* raw_virial,
                          const double fact_e,
                          const double fact_f,
                          const double fact_v,
                          double& potential,
                          ModuleBase::matrix& force,
                          ModuleBase::matrix& virial,
                          NepCudaPostprocessWorkspace& workspace);
#endif

} // namespace ModuleESolver

#endif // ESOLVER_NEP_POSTPROCESS_H
