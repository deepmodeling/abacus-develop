#ifndef SCCS_RESPONSE_H
#define SCCS_RESPONSE_H

#include "sccs_cavity.h"
#include "sccs_sa.h"
#include "source_base/vector3.h"

#include <vector>

namespace ModulePW
{
class PW_Basis;
}

namespace ModuleSccs
{
class CoulombOperator;

struct PolarizationSolverParameters;

struct PolarizationResult
{
    int iterations = 0;
    double residual_rms = 0.0;
    double residual_max = 0.0;
    bool warm_started = false;
    std::vector<double> potential; // Hartree; zero mean only for periodic cells
    std::vector<ModuleBase::Vector3<double>> gradient; // Hartree/Bohr; empty with lowpass
};

struct SccsResponse
{
    // Dielectric boundary; with the solvent-aware filling it is s_sa, while
    // dsolute_drho stays ds/dn of the local boundary s.
    std::vector<double> solute;
    std::vector<double> dsolute_drho;
    FilledCavity solvent_aware; // empty without the filling
    std::vector<double> epsilon;
    std::vector<ModuleBase::Vector3<double>> grad_log_epsilon;
    PolarizationResult polarization;
    // Hartree: reaction-energy derivative with respect to the dielectric
    // boundary (s_sa when filled); the discrete derivative with lowpass on an
    // open boundary, otherwise the continuum derivative L eps |grad v|^2/(8 pi)
    // with L = ln(eps_bulk).
    std::vector<double> boundary_potential;
    // Hartree: the same derivative with respect to the cavity density.
    std::vector<double> cavity_potential;
    std::vector<double> restart_potential; // unshifted solution for explicit warm starts
    double far_field_polarization_charge = 0.0; // open-boundary screened charge minus solute charge
};

// Chain a derivative with respect to the boundary to the cavity density,
// through the solvent-aware adjoint when the boundary is filled; probe_kernel
// is the one the response was solved with.
void boundary_to_rho(const SccsResponse& response,
                     const std::vector<double>& probe_kernel,
                     const ModulePW::PW_Basis& basis,
                     const std::vector<double>& boundary_potential,
                     std::vector<double>& density_potential);

// Periodic sqrt-CG response with chain-rule cavity derivatives, as Environ's
// CG solver with the sqrt preconditioner, without the solvent-aware filling;
// charge is ions minus electrons.
// All ranks in the PW pool must call together with the same configuration.
// initial_potential is empty on all ranks for a cold start, otherwise a local grid.
// A solver breakdown or missed tolerance stops the run with WARNING_QUIT.
// No state or INPUT globals are read or retained.
void solve_response(const std::vector<double>& cavity_density,
                    const std::vector<double>& solute_charge,
                    const CavityParameters& cavity,
                    const PolarizationSolverParameters& solver,
                    const std::vector<double>& initial_potential,
                    const ModulePW::PW_Basis& basis,
                    double tpiba,
                    SccsResponse& result);

// Use the same operator for the dielectric response and vacuum subtraction.
// cavity.solvent_aware decides on every pool rank whether the cavity is
// filled; probe_kernel is then its sa_probe_kernel, and is unused (may be
// empty) without the filling.
void solve_response(const std::vector<double>& cavity_density,
                    const std::vector<double>& solute_charge,
                    const CavityParameters& cavity,
                    const std::vector<double>& probe_kernel,
                    const PolarizationSolverParameters& solver,
                    const std::vector<double>& initial_potential,
                    const ModulePW::PW_Basis& basis,
                    double tpiba,
                    CoulombOperator& coulomb,
                    SccsResponse& result);
} // namespace ModuleSccs

#endif
