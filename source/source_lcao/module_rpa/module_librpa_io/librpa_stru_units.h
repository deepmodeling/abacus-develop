#ifndef LIBRPA_STRU_UNITS_H
#define LIBRPA_STRU_UNITS_H

namespace RpaLriDetail
{
struct LibRpaStruUnitScales
{
    double real_space_bohr;
    double reciprocal_space_bohr_inv;
};

LibRpaStruUnitScales librpa_stru_unit_scales(double lat0_bohr);
} // namespace RpaLriDetail

#endif
