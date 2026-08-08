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
    double energy_sum = 0.0;
#pragma omp parallel for reduction(+:energy_sum) schedule(static) if (nat >= 256)
    for (int i = 0; i < nat; ++i)
    {
        energy_sum += atomic_energy[i];
    }
    potential = energy_sum * fact_e;

#pragma omp parallel for schedule(static) if (nat >= 256)
    for (int i = 0; i < nat; ++i)
    {
        force(i, 0) = raw_force[i] * fact_f;
        force(i, 1) = raw_force[i + nat] * fact_f;
        force(i, 2) = raw_force[i + 2 * nat] * fact_f;
    }

    double v0 = 0.0;
    double v1 = 0.0;
    double v2 = 0.0;
    double v3 = 0.0;
    double v4 = 0.0;
    double v5 = 0.0;
    double v6 = 0.0;
    double v7 = 0.0;
    double v8 = 0.0;
#pragma omp parallel for reduction(+:v0, v1, v2, v3, v4, v5, v6, v7, v8) schedule(static) if (nat >= 256)
    for (int i = 0; i < nat; ++i)
    {
        v0 += raw_virial[i];
        v1 += raw_virial[nat + i];
        v2 += raw_virial[2 * nat + i];
        v3 += raw_virial[3 * nat + i];
        v4 += raw_virial[4 * nat + i];
        v5 += raw_virial[5 * nat + i];
        v6 += raw_virial[6 * nat + i];
        v7 += raw_virial[7 * nat + i];
        v8 += raw_virial[8 * nat + i];
    }
    const double virial_sum[9] = {v0, v1, v2, v3, v4, v5, v6, v7, v8};

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            virial(i, j) = virial_sum[3 * i + j] * fact_v;
        }
    }
}

} // namespace ModuleESolver
