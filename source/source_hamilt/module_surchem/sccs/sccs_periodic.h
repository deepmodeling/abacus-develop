#ifndef SCCS_PERIODIC_H
#define SCCS_PERIODIC_H

#include "sccs_cavity.h"
#include "sccs_poisson.h"

namespace ModulePW
{
class PW_Basis;
}

namespace ModuleSccs
{

struct PeriodicSccsResult
{
    std::vector<double> solute;
    std::vector<double> dsolute_drho;
    std::vector<double> epsilon;
    std::vector<double> depsilon_drho;
    std::vector<ModuleBase::Vector3<double>> grad_log_epsilon;
    PolarizationResult polarization;
};

PeriodicSccsResult solve_periodic_sccs(
    const std::vector<double>& cavity_density,
    const std::vector<double>& solute_charge,
    const CavityParameters& cavity_parameters,
    const PolarizationSolverParameters& solver_parameters,
    const std::vector<double>& initial_polarization_charge,
    const ModulePW::PW_Basis& basis,
    double tpiba,
    int pool_process_count);

PeriodicSccsResult solve_sccs_response(
    const std::vector<double>& cavity_density,
    const std::vector<double>& solute_charge,
    const CavityParameters& cavity_parameters,
    const PolarizationSolverParameters& solver_parameters,
    const std::vector<double>& initial_polarization_charge,
    const ModulePW::PW_Basis& basis,
    double tpiba,
    const CoulombOperator& coulomb,
    const PolarizationReduction& reduction);

} // namespace ModuleSccs

#endif
