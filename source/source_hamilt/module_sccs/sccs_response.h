#ifndef SCCS_RESPONSE_H
#define SCCS_RESPONSE_H

#include "sccs_cavity.h"
#include "source_base/vector3.h"

#include <string>
#include <vector>

namespace ModulePW
{
class PW_Basis;
}

namespace ModuleSccs
{
struct PolarizationSolverParameters;

struct PolarizationResult
{
    int iterations = 0;
    double residual_rms = 0.0;
    double residual_max = 0.0;
    bool warm_started = false;
    std::vector<double> potential; // Hartree, shifted to periodic zero mean
    std::vector<ModuleBase::Vector3<double>> gradient; // Hartree/Bohr
};

struct SccsResponse
{
    std::vector<double> solute;
    std::vector<double> dsolute_drho;
    std::vector<double> epsilon;
    std::vector<double> depsilon_drho;
    std::vector<ModuleBase::Vector3<double>> grad_log_epsilon;
    PolarizationResult polarization;
    std::vector<double> cavity_potential; // -epsilon' |grad v|^2/(8 pi), Hartree
    std::vector<double> restart_potential; // unshifted solution for explicit warm starts
};

// Original chain-derivative sqrt-CG periodic response; charge is ions minus electrons.
// All ranks in the PW pool must call together with the same configuration.
// initial_potential is empty on all ranks for a cold start, otherwise a local grid.
// Failure leaves result unchanged. No state or INPUT globals are read or retained.
bool solve_sccs_response(const std::vector<double>& cavity_density,
                         const std::vector<double>& solute_charge,
                         const CavityParameters& cavity,
                         const PolarizationSolverParameters& solver,
                         const std::vector<double>& initial_potential,
                         const ModulePW::PW_Basis& basis,
                         double tpiba,
                         SccsResponse& result,
                         std::string& error);
} // namespace ModuleSccs

#endif
