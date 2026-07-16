#ifndef MAKOV_PAYNE_H
#define MAKOV_PAYNE_H

class UnitCell;
class Charge;

namespace elecstate
{

struct MakovPayneResult
{
    double charge = 0.0;
    double first_order = 0.0;
    double second_order = 0.0;
    double total = 0.0;
    double vacuum_level = 0.0;
    bool has_vacuum_level = false;
};

/// Calculate and print the Makov-Payne isolated-cell correction for cubic cells.
/// All energies are in Ry; the optional corrected vacuum level is in eV.
MakovPayneResult makov_payne_correction(const UnitCell& ucell,
                                         const Charge& charge,
                                         const double* v_elecstat = nullptr,
                                         bool print = true);

} // namespace elecstate

#endif // MAKOV_PAYNE_H
