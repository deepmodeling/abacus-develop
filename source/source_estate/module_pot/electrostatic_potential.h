#ifndef ABACUS_ELECTROSTATIC_POTENTIAL_H
#define ABACUS_ELECTROSTATIC_POTENTIAL_H

#include <vector>
class UnitCell;
class Charge;
class surchem;
struct Input_para;
namespace ModulePW { class PW_Basis; }

namespace elecstate
{
// Same local electrostatic reference for cube output and electronic observables.
std::vector<double> electrostatic_potential(const UnitCell& cell, ModulePW::PW_Basis& basis,
                                           const Charge& charge, const double* fixed,
                                           const surchem& solvent, const Input_para& input,
                                           const std::vector<double>* pcc_potential);
double planar_potential_average(const UnitCell& cell, const ModulePW::PW_Basis& basis,
                                const std::vector<double>& potential, int axis,
                                double minimum_A, double maximum_A);
}
#endif
