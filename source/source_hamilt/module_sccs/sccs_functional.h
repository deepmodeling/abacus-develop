#ifndef SCCS_FUNCTIONAL_H
#define SCCS_FUNCTIONAL_H

#include <vector>

namespace ModulePW
{
class PW_Basis;
}

namespace ModuleSccs
{
class CoulombOperator;

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
    std::vector<double> cavity_potential; // derivative w.r.t. cavity density, Hartree
};

// Collective over the PW pool.
void eval_functional(const std::vector<double>& charge,
                      const SccsResponse& response,
                      const SccsConfig& config,
                      const ModulePW::PW_Basis& basis,
                      double tpiba,
                      FunctionalResult& result);

// Use the same operator for the dielectric response and vacuum subtraction.
void eval_functional(const std::vector<double>& charge,
                      const SccsResponse& response,
                      const SccsConfig& config,
                      const ModulePW::PW_Basis& basis,
                      double tpiba,
                      CoulombOperator& coulomb,
                      FunctionalResult& result);
} // namespace ModuleSccs

#endif
