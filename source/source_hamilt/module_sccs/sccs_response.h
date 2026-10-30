#ifndef SCCS_RESPONSE_H
#define SCCS_RESPONSE_H

#include "sccs_cavity.h"
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
    std::vector<double> solute;
    std::vector<double> dsolute_drho;
    std::vector<double> epsilon;
    std::vector<double> depsilon_drho;
    std::vector<ModuleBase::Vector3<double>> grad_log_epsilon;
    PolarizationResult polarization;
    // Hartree: discrete cavity derivative with lowpass, otherwise the continuum
    // derivative -epsilon' |grad v|^2/(8 pi).
    std::vector<double> cavity_potential;
    std::vector<double> restart_potential; // unshifted solution for explicit warm starts
    double far_field_polarization_charge = 0.0; // open-boundary screened charge minus solute charge
};

// Periodic sqrt-CG response with chain-rule cavity derivatives, as Environ's
// CG solver with the sqrt preconditioner; charge is ions minus electrons.
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
void solve_response(const std::vector<double>& cavity_density,
                    const std::vector<double>& solute_charge,
                    const CavityParameters& cavity,
                    const PolarizationSolverParameters& solver,
                    const std::vector<double>& initial_potential,
                    const ModulePW::PW_Basis& basis,
                    double tpiba,
                    CoulombOperator& coulomb,
                    SccsResponse& result);
} // namespace ModuleSccs

#endif
