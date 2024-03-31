#include "module_base/vector3.h"
/// index structure for grid integral module
/// in ABACUS, this index is stored for tracing:
/// 1. starting row and column index (mu, nu)
/// 2. R distance from atom 1 and atom2 (dR)
/// 3. number of orbitals for atom1 and atom2 (nw1, nw2)
namespace gridIntegral
{

struct gridIndex
{
    int nnrg;
    int mu;
    int nu;
    ModuleBase::Vector3<int> dR;
    int nw1;
    int nw2;
};

} // namespace gridIntegral