#ifndef PDOS_ORBITALS_H
#define PDOS_ORBITALS_H

#include <vector>

class Atom;
namespace ModuleBase
{
class matrix;
}

namespace ModuleIO
{
struct PdosOrbital
{
    int l;
    int m; // physical magnetic quantum number
    std::vector<int> zeta_orbitals;

    double value(const ModuleBase::matrix& pdos, int first, int energy, bool spinor) const;
};

std::vector<PdosOrbital> pdos_orbitals(const Atom& atom);
}
#endif
