#ifndef ABACUS_ELPOT_H
#define ABACUS_ELPOT_H

#include <utility>
#include <vector>
class UnitCell;
class Charge;
class surchem;
struct Input_para;
namespace ModulePW
{
class PW_Basis;
}

namespace elecstate
{
class Potential;
struct PotentialWindow;
double electrostatic_reference(
    const UnitCell &cell,
    ModulePW::PW_Basis &basis,
    const Charge &charge,
    const Potential &pot,
    const surchem &solvent,
    const Input_para &input,
    const PotentialWindow &window);

// Same local electrostatic reference for cube output and electronic observables.
std::vector<double> electrostatic_potential(
    const UnitCell &cell,
    ModulePW::PW_Basis &basis,
    const Charge &charge,
    const double *fixed,
    const surchem &solvent,
    const Input_para &input,
    const std::vector<double> *pcc_potential,
    const std::vector<double> *solvent_potential);
// Validate cell bounds and return the normal cell height in angstrom.
double checked_window_height(
    const UnitCell &cell,
    int axis,
    double minimum_A,
    double maximum_A);
// Return the nonempty half-open interval of global FFT plane indices.
std::pair<int, int> potential_window_planes(
    const UnitCell &cell,
    const ModulePW::PW_Basis &basis,
    int axis,
    double minimum_A,
    double maximum_A);
double planar_potential_average(
    const UnitCell &cell,
    const ModulePW::PW_Basis &basis,
    const std::vector<double> &potential,
    int axis,
    double minimum_A,
    double maximum_A);
} // namespace elecstate
#endif
