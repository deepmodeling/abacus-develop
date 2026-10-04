#include "charge_moments.h"

#include <cmath>

namespace elecstate
{

bool charge_moments(const double* charges,
                     const ModuleBase::Vector3<double>* relative_positions,
                     const int count,
                     const double weight,
                     ChargeMoments& moments,
                     std::string& error)
{
    error.clear();
    if (count < 0 || (count > 0 && (charges == nullptr || relative_positions == nullptr))
        || !std::isfinite(weight) || weight <= 0.0)
    {
        error = "charge moments require valid sample storage and a positive finite weight";
        return false;
    }
    ChargeMoments candidate;
    for (int index = 0; index < count; ++index)
    {
        const double charge = charges[index] * weight;
        const ModuleBase::Vector3<double>& relative = relative_positions[index];
        if (!std::isfinite(charge) || !std::isfinite(relative.x)
            || !std::isfinite(relative.y) || !std::isfinite(relative.z))
        {
            error = "charge moments require finite charges and relative positions";
            return false;
        }
        candidate.charge += charge;
        candidate.dipole += relative * charge;
        candidate.second_moment += charge * relative.norm2();
    }
    if (!std::isfinite(candidate.charge) || !std::isfinite(candidate.dipole.x)
        || !std::isfinite(candidate.dipole.y) || !std::isfinite(candidate.dipole.z)
        || !std::isfinite(candidate.second_moment))
    {
        error = "charge moment accumulation must be finite";
        return false;
    }
    moments = candidate;
    return true;
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
