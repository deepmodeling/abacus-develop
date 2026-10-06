#ifndef SCCS_FUNCTIONAL_H
#define SCCS_FUNCTIONAL_H

#include <vector>

namespace ModulePW
{
class PW_Basis;
}

namespace ModuleSccs
{
struct SccsConfig;
struct SccsResponse;

struct FunctionalResult
{
    double reaction_energy = 0.0; // Hartree
    double surface_energy = 0.0;
    double volume_energy = 0.0;
    double surface = 0.0; // Bohr^2
    double volume = 0.0; // Bohr^3
    std::vector<double> reaction_potential; // dielectric minus vacuum, Hartree
    std::vector<double> electron_potential; // derivative of all three energy terms
};

// Collective over the PW pool; a non-finite energy stops the run.
void evaluate_functional(const std::vector<double>& charge,
                          const SccsResponse& response,
                          const SccsConfig& config,
                          const ModulePW::PW_Basis& basis,
                          double tpiba,
                          FunctionalResult& result);
} // namespace ModuleSccs

#endif
