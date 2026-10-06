#include "charge_moments.h"

namespace elecstate
{

ChargeMoments charge_moments(const double* charges,
                             const ModuleBase::Vector3<double>* relative_positions,
                             const int count,
                             const double weight)
{
    ChargeMoments moments;
    for (int index = 0; index < count; ++index)
    {
        const double charge = charges[index] * weight;
        const ModuleBase::Vector3<double>& relative = relative_positions[index];
        moments.charge += charge;
        moments.dipole += relative * charge;
        moments.second_moment += charge * relative.norm2();
    }
    return moments;
}

ChargeMoments add_charge_moments(const ChargeMoments& left, const ChargeMoments& right)
{
    ChargeMoments result;
    result.charge = left.charge + right.charge;
    result.dipole = left.dipole + right.dipole;
    result.second_moment = left.second_moment + right.second_moment;
    return result;
}

} // namespace elecstate
