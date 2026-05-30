#include "operator_force_stress_utils.h"
#include "source_base/parallel_reduce.h"

namespace OperatorForceStress {

void finalize_force_stress(
    bool cal_force,
    bool cal_stress,
    const UnitCell* ucell,
    const std::vector<double>& stress_tmp,
    ModuleBase::matrix& force,
    ModuleBase::matrix& stress,
    double force_factor,
    double stress_factor)
{
    if (cal_force)
    {
#ifdef __MPI
        Parallel_Reduce::reduce_all(force.c, force.nr * force.nc);
#endif
        // Apply factor of 2 for Hermitian matrix
        for (int i = 0; i < force.nr * force.nc; i++)
        {
            force.c[i] *= force_factor;
        }
    }

    if (cal_stress)
    {
#ifdef __MPI
        Parallel_Reduce::reduce_all(const_cast<double*>(stress_tmp.data()), 6);
#endif
        const double weight = ucell->lat0 / ucell->omega;
        for (int i = 0; i < 6; i++)
        {
            stress.c[i] = stress_tmp[i] * weight * stress_factor;
        }
        // Rearrange to 3x3 matrix format
        rearrange_stress_matrix(stress);
    }
}

} // namespace OperatorForceStress
