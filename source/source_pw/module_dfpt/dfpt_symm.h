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

/// One (k, k+q) response pair representative with its orbit weight.
struct KQPair
{
    ModuleBase::Vector3<double> k_frac;  ///< representative mesh k, wrapped to [0,1)
    ModuleBase::Vector3<double> kq_frac; ///< k + q, wrapped to [0,1)
    int orbit_size;                      ///< number of mesh points this pair represents
};

/**
 * @brief Representatives of the (k, k+q) response pairs over a full mesh (B1).
 *
 * Reduces the full k mesh to orbit representatives of the UNITARY little
 * co-group of q (the s = +1 branch: QE's nsymq), and returns one pair
 * (k, k + q) per orbit with the orbit size as weight. The plain weighted
 * sum over these representatives is the input the G_q projection of the
 * later symmetrization step closes on (P(A) = F requires the list to be
 * reduced with respect to the same group the projector runs over).
 *
 * The orbits are generated in cartesian k space through the rotation
 * images R = latvec^-1 G latvec, which makes the orbit PARTITION immune
 * to the forward/backward gmatrix conventions ({R} and {R^T} generate
 * the same orbits because the cartesian images are orthogonal); the
 * Gamma case reproduces the standard k-mesh IBZ reduction as the
 * convention anchor.
 *
 * Inert infrastructure: nothing consumes it yet.
 */
std::vector<KQPair> kq_pair_reps(
    const ModuleBase::Matrix3 &latvec,
    const ModuleBase::Matrix3 *gmatrix,
    int nops,
    const ModuleBase::Vector3<double> &q_frac,
    const std::vector<ModuleBase::Vector3<double>> &mesh_frac);

} // namespace ModuleDFPT

#endif // DFPT_SYMM_H
