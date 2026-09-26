#ifndef SOURCE_DFTB_SK_BLOCK_H
#define SOURCE_DFTB_SK_BLOCK_H

#include "source_dftb/skf_data.h"

#include <array>

namespace ModuleDFTB
{

/**
 * Legacy v1.0 SKF column indices (zero-based, as stored in the file).
 * The file order is dd(sigma,pi,delta), pd(sigma,pi), pp(sigma,pi),
 * sd(sigma), sp(sigma), ss(sigma).
 */
enum class SkfChannel : std::size_t
{
    dd_sigma = 0,
    dd_pi = 1,
    dd_delta = 2,
    pd_sigma = 3,
    pd_pi = 4,
    pp_sigma = 5,
    pp_pi = 6,
    sd_sigma = 7,
    sp_sigma = 8,
    ss_sigma = 9
};

/**
 * One real 4x4 diatomic H/S block for one s and one p shell on each atom.
 * Rows are orbitals on atom B; columns are orbitals on atom A. Orbital order
 * is s, p_y, p_z, p_x, matching DFTB+'s real-harmonic Slater-Koster code.
 * Cartesian derivatives are with respect to R = r_B - r_A, in atomic units.
 */
struct SkfSpBlockEvaluation
{
    std::array<double, 16> hamiltonian{};
    std::array<double, 16> overlap{};
    std::array<std::array<double, 16>, 3> d_hamiltonian_dR{};
    std::array<std::array<double, 16>, 3> d_overlap_dR{};
};

/**
 * Build an sp-only pair block from the directed A-B and B-A legacy SKF files.
 * Each pair file must contain matching grid metadata; its H/S tables are
 * allowed (and expected) to differ by direction.
 */
SkfSpBlockEvaluation evaluate_sp_pair_block(const SkfData& ab,
                                             const SkfData& ba,
                                             const std::array<double, 3>& displacement_bohr);

} // namespace ModuleDFTB

#endif // SOURCE_DFTB_SK_BLOCK_H
