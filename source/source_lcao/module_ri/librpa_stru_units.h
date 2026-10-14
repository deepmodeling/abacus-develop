#ifndef LIBRPA_STRU_UNITS_H
#define LIBRPA_STRU_UNITS_H

#include "source_base/constants.h"

#include <cmath>
#include <stdexcept>

namespace RpaLriDetail
{
struct LibRpaStruUnitScales
{
    double real_space_bohr;
    double reciprocal_space_bohr_inv;
};

inline LibRpaStruUnitScales librpa_stru_unit_scales(const double lat0_bohr)
{
    if (!std::isfinite(lat0_bohr) || lat0_bohr <= 0.0)
    {
        throw std::invalid_argument("LibRPA stru_out requires a positive finite lattice constant.");
    }
    return {lat0_bohr, ModuleBase::TWO_PI / lat0_bohr};
}
} // namespace RpaLriDetail

#endif
