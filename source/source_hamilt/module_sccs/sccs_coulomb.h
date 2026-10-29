#ifndef SCCS_COULOMB_H
#define SCCS_COULOMB_H

#include <vector>

namespace ModuleSccs
{
// Potential in Hartree for signed charge density in e/Bohr^3. Applications
// are collective over the PW pool.
class CoulombOperator
{
public:
    virtual ~CoulombOperator() {}
    virtual bool has_correction() const = 0;
    virtual void apply_potential(const std::vector<double>& charge, std::vector<double>& potential) = 0;
};
} // namespace ModuleSccs

#endif
