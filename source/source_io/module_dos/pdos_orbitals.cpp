#include "pdos_orbitals.h"
#include "source_base/matrix.h"
#include "source_cell/atom_spec.h"

#include <cmath>
#include <utility>

std::vector<ModuleIO::PdosOrbital> ModuleIO::pdos_orbitals(const Atom& atom)
{
    std::vector<PdosOrbital> orbitals;
    for (int l = 0; l <= atom.nwl; ++l)
    {
        for (int im = 0; im <= 2 * l; ++im)
        {
            const int m = im % 2 == 0 ? -im / 2 : (im + 1) / 2;
            PdosOrbital orbital{l, m, {}};
            for (int j = 0; j < atom.nw; ++j)
            {
                if (atom.iw2l[j] == l && atom.iw2m[j] == im)
                {
                    orbital.zeta_orbitals.push_back(j);
                }
            }
            if (!orbital.zeta_orbitals.empty())
            {
                orbitals.push_back(std::move(orbital));
            }
        }
    }
    return orbitals;
}

double ModuleIO::PdosOrbital::value(const ModuleBase::matrix& pdos, int first, int energy, bool spinor) const
{
    double sum = 0.0;
    for (int j : zeta_orbitals)
    {
        if (spinor)
        {
            const int row = first + 2 * j;
            sum += pdos(row, energy) + pdos(row + 1, energy);
        }
        else
        {
            sum += pdos(first + j, energy);
        }
    }
    return std::abs(sum) < 1.0e-6 ? 0.0 : sum;
}
