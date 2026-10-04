#ifndef SCCS_PW_COULOMB_H
#define SCCS_PW_COULOMB_H

#include <complex>
#include <string>
#include <vector>

namespace ModulePW
{
class PW_Basis;
}

namespace ModuleSccs
{
// All calls are collective over the ABACUS pool used to initialize basis.
bool validate_pw_grid(const ModulePW::PW_Basis& basis, double tpiba, std::string& error);
bool validate_grid_values(const std::vector<double>& values,
                          const ModulePW::PW_Basis& basis,
                          std::string& error);

// Potential in Hartree for a signed charge density in e/Bohr^3.
// G=0 is removed: this is the periodic Coulomb inverse with a uniform background.
class PeriodicCoulombOperator
{
public:
    PeriodicCoulombOperator(const ModulePW::PW_Basis& basis, double tpiba);
    bool apply_potential(const std::vector<double>& charge,
                         std::vector<double>& potential,
                         std::string& error);

private:
    const ModulePW::PW_Basis& basis_;
    const double tpiba_;
    // Scratch only, overwritten on each call; PW_Basis FFTs are also non-reentrant.
    std::vector<std::complex<double>> reciprocal_work_;
};
} // namespace ModuleSccs

#endif
