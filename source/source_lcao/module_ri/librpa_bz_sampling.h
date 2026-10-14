#ifndef LIBRPA_BZ_SAMPLING_H
#define LIBRPA_BZ_SAMPLING_H

#include "source_base/constants.h"
#include "source_base/vector3.h"

#include <cmath>
#include <stdexcept>

namespace RpaLriDetail
{
inline ModuleBase::Vector3<double> librpa_cartesian_kvec(const ModuleBase::Vector3<double>& kvec_c,
                                                         const double lat0_bohr)
{
    if (!std::isfinite(lat0_bohr) || lat0_bohr <= 0.0)
    {
        throw std::invalid_argument("LibRPA BZ output requires a positive lattice constant.");
    }
    return kvec_c * (ModuleBase::TWO_PI / lat0_bohr);
}
} // namespace RpaLriDetail

#endif
