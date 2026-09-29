#ifndef SCCS_RESPONSE_H
#define SCCS_RESPONSE_H

#include "sccs_cavity.h"
#include "sccs_poisson.h"

namespace ModulePW
{
class PW_Basis;
}

namespace ModuleSurchem
{
class ChargeReduction;
}

namespace ModuleSccs
{

// Dielectric response of one SCCS evaluation: the cavity fields and the
// sqrt-CG solution for the solute charge.
struct SccsResponse
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

// ENVIRON dielectric_of_potential polarization density,
// grad(ln eps).grad(v)/(4 pi) + q (1/eps - 1). On a finite grid it need not
// equal -laplacian(v)/(4 pi) - q; ABACUS uses it only for PCC diagnostics.
std::vector<double> continuum_polarization_charge(
    const std::vector<double>& solute_charge,
    const SccsResponse& response);

// ENVIRON sqrt-preconditioned CG for every boundary: the preconditioner uses
// coulomb (periodic or PCC-corrected). It stops on the RMS and maximum charge
// residual and warm-starts from initial_potential (previous solution, or empty)
// when that helps.
SccsResponse solve_sccs_response(
    const std::vector<double>& cavity_density,
    const std::vector<double>& solute_charge,
    const CavityParameters& cavity_parameters,
    const PolarizationSolverParameters& solver_parameters,
    const std::vector<double>& initial_potential,
    const ModulePW::PW_Basis& basis,
    double tpiba,
    const CoulombOperator& coulomb,
    const ModuleSurchem::ChargeReduction& reduction);

} // namespace ModuleSccs

#endif
