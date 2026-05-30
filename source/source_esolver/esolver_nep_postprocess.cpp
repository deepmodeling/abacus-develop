#include "esolver_nep_postprocess.h"

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
                         ModuleBase::matrix& virial)
{
    potential = 0.0;
    for (int i = 0; i < nat; ++i)
    {
        potential += atomic_energy[i] * fact_e;
        force(i, 0) = raw_force[i] * fact_f;
        force(i, 1) = raw_force[i + nat] * fact_f;
        force(i, 2) = raw_force[i + 2 * nat] * fact_f;
    }

    double virial_sum[9] = {0.0};
    for (int j = 0; j < 9; ++j)
    {
        const int offset = j * nat;
        for (int i = 0; i < nat; ++i)
        {
            virial_sum[j] += raw_virial[offset + i];
        }
    }

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            virial(i, j) = virial_sum[3 * i + j] * fact_v;
        }
    }
}

} // namespace ModuleESolver
