#ifndef LIBRPA_BZ_SAMPLING_H
#define LIBRPA_BZ_SAMPLING_H

#include "source_base/constants.h"
#include "source_base/vector3.h"

#include <cmath>
#include <stdexcept>

namespace RpaLriDetail
{
struct LibRpaStoredQIndex
{
    int coulomb_irreducible_index;
    int representative_scf_index;
};

inline ModuleBase::Vector3<double> librpa_cartesian_kvec(const ModuleBase::Vector3<double>& kvec_c,
                                                         const double lat0_bohr)
{
    if (!std::isfinite(lat0_bohr) || lat0_bohr <= 0.0)
    {
        throw std::invalid_argument("LibRPA BZ output requires a positive lattice constant.");
    }
    return kvec_c * (ModuleBase::TWO_PI / lat0_bohr);
}

inline int librpa_stored_coulomb_q_count(const int stored_kpoint_count)
{
    if (stored_kpoint_count <= 0)
    {
        throw std::invalid_argument("LibRPA BZ output requires stored k points.");
    }
    return stored_kpoint_count;
}

inline LibRpaStoredQIndex librpa_stored_q_index(const int stored_kpoint_index,
                                                const int stored_kpoint_count)
{
    if (stored_kpoint_index < 0 || stored_kpoint_index >= stored_kpoint_count)
    {
        throw std::invalid_argument("LibRPA BZ output found an invalid stored k-point index.");
    }
    const int one_based_index = stored_kpoint_index + 1;
    return {one_based_index, one_based_index};
}
} // namespace RpaLriDetail

#endif
