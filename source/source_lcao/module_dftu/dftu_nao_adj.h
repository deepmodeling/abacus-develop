#ifndef DFTU_NAO_ADJ_H
#define DFTU_NAO_ADJ_H

/// @file dftu_nao_adj.h
/// @brief Structure-dependent precomputation for the DFT+U LCAO operator:
///        Hubbard-atom adjacent lists and the <phi|chi_m> overlap table.

#include <unordered_map>
#include <vector>

#include "dftu_nao_ijr.h" // DFTU_LCAO::NlmTot

class UnitCell;
class Plus_U_Base;
class Grid_Driver;
class AdjacentAtomInfo;
class Parallel_Orbitals;
class TwoCenterIntegrator;

namespace DFTU_LCAO
{

/**
 * @brief build the adjacent-atom lists for all Hubbard atoms.
 *
 * For every atom whose type carries a Hubbard-U channel, find its neighbor
 * atoms within orb_cutoff + onsite_radius and filter the AdjacentAtomInfo
 * accordingly.
 *
 * @param ucell         [in] unit cell
 * @param dftu          [in] DFT+U base object (per-type U channels)
 * @param gridD         [in] grid driver for neighbor search
 * @param orb_cutoff    [in] orbital cutoff radius per atom type
 * @param onsite_radius [in] onsite projector cutoff radius
 * @return one AdjacentAtomInfo per Hubbard atom, in atom order
 */
std::vector<AdjacentAtomInfo> build_adjacent_atoms(const UnitCell* ucell,
                                                   Plus_U_Base* dftu,
                                                   const Grid_Driver* gridD,
                                                   const std::vector<double>& orb_cutoff,
                                                   const double onsite_radius);

/**
 * @brief calculate the <phi|alpha^I> overlap values for all Hubbard atoms.
 *
 * The result is indexed by global atom index (empty entries for non-Hubbard
 * atoms) and reused in the HR/occupation accumulation.
 *
 * @param ucell    [in] unit cell
 * @param dftu     [in] DFT+U base object (per-type U channels)
 * @param intor    [in] two-center integrator for <phi|chi_m>
 * @param adjs_all [in] adjacent atom info for all Hubbard atoms
 * @param pv       [in] parallel-orbitals descriptor providing local index maps
 * @return the overlap table, one outer entry per atom (ucell->nat entries)
 */
NlmTot cal_nlm_all(const UnitCell& ucell,
                   const Plus_U_Base& dftu,
                   const TwoCenterIntegrator& intor,
                   const std::vector<AdjacentAtomInfo>& adjs_all,
                   const Parallel_Orbitals& pv);

} // namespace DFTU_LCAO

#endif // DFTU_NAO_ADJ_H
