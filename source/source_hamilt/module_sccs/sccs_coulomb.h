#ifndef SCCS_COULOMB_H
#define SCCS_COULOMB_H

#include <vector>

namespace ModuleSccs
{
struct CoulombTransformCounts
{
    int forward_calls = 0;
    int inverse_calls = 0;
};

// Potential in Hartree for signed charge density in e/Bohr^3. Applications
// are collective over the PW pool.
class CoulombOperator
{
public:
    virtual ~CoulombOperator() {}
    virtual bool has_correction() const = 0;
    virtual CoulombTransformCounts fft_counts() const { return {}; }
    virtual void apply_potential(const std::vector<double>& charge, std::vector<double>& potential) = 0;
};
} // namespace ModuleSccs

#endif
