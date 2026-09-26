#include "sccs_periodic.h"

#include "sccs_pw_coulomb.h"
#include "sccs_pw_reduction.h"

#include "source_basis/module_pw/pw_basis.h"

#include <cmath>
#include <stdexcept>

namespace ModuleSccs
{

PeriodicSccsResult solve_periodic_sccs(
    const std::vector<double>& cavity_density,
    const std::vector<double>& solute_charge,
    const CavityParameters& cavity_parameters,
    const PolarizationSolverParameters& solver_parameters,
    const std::vector<double>& initial_polarization_charge,
    const ModulePW::PW_Basis& basis,
    const double tpiba,
    const int pool_process_count)
{
    const PeriodicCoulombOperator coulomb(basis, tpiba);
    const PoolPolarizationReduction reduction(pool_process_count);
    return solve_sccs_response(cavity_density,
                               solute_charge,
                               cavity_parameters,
                               solver_parameters,
                               initial_polarization_charge,
                               basis,
                               tpiba,
                               coulomb,
                               reduction);
}

PeriodicSccsResult solve_sccs_response(
    const std::vector<double>& cavity_density,
    const std::vector<double>& solute_charge,
    const CavityParameters& cavity_parameters,
    const PolarizationSolverParameters& solver_parameters,
    const std::vector<double>& initial_polarization_charge,
    const ModulePW::PW_Basis& basis,
    const double tpiba,
    const CoulombOperator& coulomb,
    const PolarizationReduction& reduction)
{
    if (cavity_density.size() != static_cast<std::size_t>(basis.nrxx)
        || solute_charge.size() != cavity_density.size())
    {
        throw std::invalid_argument("SCCS arrays must match the local PW real-space grid");
    }
    validate_cavity_parameters(cavity_parameters);
    validate_polarization_solver_parameters(solver_parameters);

    PeriodicSccsResult result;
    result.solute.resize(cavity_density.size());
    result.dsolute_drho.resize(cavity_density.size());
    result.epsilon.resize(cavity_density.size());
    result.depsilon_drho.resize(cavity_density.size());
    std::vector<double> log_epsilon(cavity_density.size());
    for (std::size_t index = 0; index < cavity_density.size(); ++index)
    {
        const CavityPoint point = evaluate_cavity(cavity_density[index], cavity_parameters);
        result.solute[index] = point.solute;
        result.dsolute_drho[index] = point.dsolute_drho;
        result.epsilon[index] = point.epsilon;
        result.depsilon_drho[index] = point.depsilon_drho;
        log_epsilon[index] = std::log(point.epsilon);
    }

    result.grad_log_epsilon = periodic_gradient(log_epsilon, basis, tpiba);
    result.polarization = solve_polarization(solute_charge,
                                             result.epsilon,
                                             result.grad_log_epsilon,
                                             initial_polarization_charge,
                                             solver_parameters,
                                             coulomb,
                                             reduction);
    return result;
}

} // namespace ModuleSccs
