#include "sccs_moments.h"

#include "source_base/parallel_reduce.h"
#include "source_basis/module_pw/pw_basis.h"

namespace ModuleSccs
{
elecstate::ChargeMoments pool_moments(const std::vector<double>& charge,
                                      const std::vector<ModuleBase::Vector3<double>>& positions,
                                      const ModulePW::PW_Basis& basis)
{
    const double dv = basis.omega / basis.nxyz;
    const int count = charge.size();
    const double* charge_data = charge.data();
    const ModuleBase::Vector3<double>* position_data = positions.data();
    elecstate::ChargeMoments moments = elecstate::charge_moments(charge_data, position_data, count, dv);
    double values[5] = {moments.charge, moments.dipole.x, moments.dipole.y, moments.dipole.z,
                        moments.second_moment};
    Parallel_Reduce::reduce_pool(values, 5);
    moments.charge = values[0];
    moments.dipole = ModuleBase::Vector3<double>(values[1], values[2], values[3]);
    moments.second_moment = values[4];
    return moments;
}
} // namespace ModuleSccs
