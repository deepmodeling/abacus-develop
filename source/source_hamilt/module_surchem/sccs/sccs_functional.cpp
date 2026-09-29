#include "sccs_functional.h"

#include <cmath>
#include <stdexcept>

namespace ModuleSccs
{

ElectrostaticFunctionalResult evaluate_electrostatic_functional(
    const std::vector<double>& solute_charge,
    const ElectrostaticField& dielectric_field,
    const ElectrostaticField& vacuum_field,
    const std::vector<double>& cavity_potential,
    const double volume_element,
    const ChargeReduction& reduction)
{
    const std::size_t size = solute_charge.size();
    if (size == 0 || dielectric_field.potential.size() != size
        || vacuum_field.potential.size() != size || cavity_potential.size() != size)
    {
        throw std::invalid_argument("SCCS electrostatic functional arrays must have the same non-zero size");
    }
    if (!std::isfinite(volume_element) || volume_element <= 0.0)
    {
        throw std::invalid_argument("SCCS electrostatic functional requires a positive finite volume element");
    }

    ElectrostaticFunctionalResult result;
    result.reaction_potential.resize(size);
    result.electron_potential.resize(size);
    for (std::size_t index = 0; index < size; ++index)
    {
        if (!std::isfinite(solute_charge[index])
            || !std::isfinite(dielectric_field.potential[index])
            || !std::isfinite(vacuum_field.potential[index])
            || !std::isfinite(cavity_potential[index]))
        {
            throw std::domain_error("SCCS electrostatic functional inputs must be finite");
        }
        const double reaction_potential
            = dielectric_field.potential[index] - vacuum_field.potential[index];
        result.reaction_potential[index] = reaction_potential;
        result.reaction_energy
            += 0.5 * solute_charge[index] * reaction_potential * volume_element;
        result.electron_potential[index] = -reaction_potential + cavity_potential[index];
    }
    result.charge_potential = result.reaction_potential;
    reduction.reduce_sum(result.reaction_energy);
    if (!std::isfinite(result.reaction_energy))
    {
        throw std::domain_error("SCCS reduced reaction energy must be finite");
    }
    return result;
}

} // namespace ModuleSccs
