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
    // Unscaled Dn from this response's cavity density. Keep it for the
    // discrete derivative, including points where epsilon' is zero.
    std::vector<ModuleBase::Vector3<double>> density_gradient;
    std::vector<ModuleBase::Vector3<double>> grad_log_epsilon;
    PolarizationResult polarization;
    // Derivative of the reaction energy with respect to the cavity density
    // through epsilon, in Ha. With the switching lowpass (PCC only) it is the
    // exact derivative of the discrete sqrt-CG energy; otherwise it is the
    // continuum -eps'|grad v|^2/(8 pi) of Environ.
    std::vector<double> cavity_potential;
    // Unshifted sqrt-CG solution: the fixed point used for the next warm start.
    std::vector<double> restart_potential;
    // Open-boundary (PCC) solutions only: the net polarization charge seen by
    // the far field, which screens the solute to int(s)/sqrt(eps_bulk) because
    // sqrt(eps) v = C_PCC(s) with s = (q - f v)/sqrt(eps). Zero when periodic.
    double far_field_polarization_charge = 0.0;
};

// Continuous dielectric source used by Environ's ionic-force path.
// On a finite grid it need not equal -laplacian(phi)/(4*pi)-q.
std::vector<double> continuum_polarization_charge(
    const std::vector<double>& solute_charge,
    const PeriodicSccsResult& response);

// initial_potential: previous sqrt-CG potential for the warm start, or empty.
PeriodicSccsResult solve_periodic_sccs(
    const std::vector<double>& cavity_density,
    const std::vector<double>& solute_charge,
    const CavityParameters& cavity_parameters,
    const PolarizationSolverParameters& solver_parameters,
    const std::vector<double>& initial_potential,
    const ModulePW::PW_Basis& basis,
    double tpiba,
    int pool_process_count);

// ENVIRON sqrt-preconditioned CG for every boundary: the preconditioner uses
// coulomb (periodic or PCC-corrected). It stops on the RMS and maximum charge
// residual and warm-starts from initial_potential (previous solution) when that helps.
PeriodicSccsResult solve_chain_sccs_response(
    const std::vector<double>& cavity_density,
    const std::vector<double>& solute_charge,
    const CavityParameters& cavity_parameters,
    const PolarizationSolverParameters& solver_parameters,
    const std::vector<double>& initial_potential,
    const ModulePW::PW_Basis& basis,
    double tpiba,
    const CoulombOperator& coulomb,
    const PolarizationReduction& reduction);

} // namespace ModuleSccs

#endif
