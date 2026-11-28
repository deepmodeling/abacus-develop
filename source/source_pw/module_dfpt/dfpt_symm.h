#ifndef DFPT_SYMM_H
#define DFPT_SYMM_H

#include "source_base/matrix3.h"
#include "source_base/vector3.h"

#include <vector>

namespace ModuleDFPT
{

/**
 * @brief Little group of a wave vector including time reversal (B1).
 *
 * Returns the indices of the space-group operations whose rotation maps
 * q onto +q or -q modulo a reciprocal-lattice vector:
 * R q = s q + G, s = +-1: the union of the unitary little co-group
 * (s = +1, QE's nsymq) and the antiunitary-eligible set (s = -1, the
 * rotations that preserve q when composed with time reversal; the
 * response at -q is the conjugate response at q). For a
 * time-reversal-invariant q the two branches coincide; for a generic q
 * the antiunitary branch still contributes inversion (theta*i fixes
 * every q), so the result always contains the identity and inversion.
 * Callers needing only the unitary co-group (the symdvscf rotation
 * average) must re-filter on the s = +1 branch.
 *
 * The reciprocal test is the dot-product integrality
 * (R q - s q) . a_i / 2pi in Z for every direct primitive vector a_i
 * (rows of latvec). This form is immune to the all-same-parity
 * subtlety of non-primitive lattices that a cartesian-integrality
 * check misses: at the fcc L point the cartesian difference
 * (1, 0, -1)*2pi/a is integral but not a reciprocal-lattice vector,
 * which would wrongly admit C2y into the little group.
 *
 * Membership is convention-immune: an operation and its inverse either
 * both fix q mod G or both do not, so the forward/backward gmatrix
 * conventions give the same index set.
 */
std::vector<int> little_group_of_q(
    const ModuleBase::Matrix3 &latvec,
    const ModuleBase::Matrix3 *gmatrix,
    int nops,
    const ModuleBase::Vector3<double> &q_cart);

} // namespace ModuleDFPT

#endif // DFPT_SYMM_H
