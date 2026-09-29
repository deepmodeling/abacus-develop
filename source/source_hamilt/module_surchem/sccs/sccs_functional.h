#ifndef SCCS_FUNCTIONAL_H
#define SCCS_FUNCTIONAL_H

#include <vector>

namespace ModuleSurchem
{
class ChargeReduction;
}

namespace ModuleSccs
{

struct ElectrostaticFunctionalResult
{
    double reaction_energy = 0.0;
    std::vector<double> reaction_potential;
    std::vector<double> electron_potential;
};

// Reaction energy 1/2 integral q (phi_dielectric - phi_vacuum), in Ha.
// The returned electronic potential is -(phi_dielectric - phi_vacuum) plus
// cavity_potential, the solver's derivative through the dielectric cavity.
ElectrostaticFunctionalResult evaluate_electrostatic_functional(
    const std::vector<double>& solute_charge,
    const std::vector<double>& dielectric_potential,
    const std::vector<double>& vacuum_potential,
    const std::vector<double>& cavity_potential,
    double volume_element,
    const ModuleSurchem::ChargeReduction& reduction);

} // namespace ModuleSccs

#endif
