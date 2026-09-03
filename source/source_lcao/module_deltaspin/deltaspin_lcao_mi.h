/**
 * @file deltaspin_lcao_mi.h
 * @brief LCAO-specific magnetic-moment calculation for DeltaSpin,
 *        expressed as free functions over ScState.
 *
 * @par Purpose
 * Holds the LCAO-only half of the DeltaSpin Mi pipeline that used to live
 * as SpinConstrain member functions:
 *   - cal_mi_lcao(): moments via the DeltaSpin operator on the density matrix
 *   - convert()/calculate_MW(): moments via the orbital multiplication matrix
 *     (alternative/debug path)
 *   - collect_MW(): accumulate mu*density-matrix contributions (ScaLAPACK)
 *
 * Keeping these as free functions decouples the LCAO path from the
 * SpinConstrain singleton: the only inputs are the constraint state
 * (ScState), the LCAO operator, the density matrix and the parallel
 * orbitals mapping.
 */
#ifndef DELTASPIN_LCAO_MI_H
#define DELTASPIN_LCAO_MI_H

#ifdef __LCAO

#include <vector>

#include "source_base/complexmatrix.h"
#include "source_base/matrix.h"
#include "source_hamilt/operator.h"

#include "deltaspin_state.h"

class Parallel_Orbitals;
namespace elecstate
{
template <typename TK, typename TR>
class DensityMatrix;
}

namespace spinconstrain
{
namespace lcao
{

/**
 * @brief Calculate atomic magnetic moments from the density matrix (LCAO).
 *
 * @details Uses the DeltaSpin operator to compute Tr(rho * mu) per atom.
 * For nspin=2, extracts only the z-component. For nspin=4, extracts
 * all three components from the interleaved 4-component spinor density matrix.
 * Results are stored in state.Mi_ (indexed by global atom index iat).
 *
 * @param state  Constraint state (Mi_ written, indexing maps read)
 * @param p_operator Base pointer to DeltaSpin<OperatorLCAO<...>>; nullptr aborts
 * @param dm     Density matrix (rho in orbital basis)
 * @param step   Current SCF iteration number (for logging)
 * @param print  Whether to print moments to ofs_running
 */
void cal_mi_lcao(ScState& state,
                 hamilt::Operator<std::complex<double>>* p_operator,
                 elecstate::DensityMatrix<std::complex<double>, double>* dm,
                 const int& step,
                 bool print = false);

/**
 * @brief Convert flat orbital matrix to nested vector [nspin][iat][iw].
 *
 * @param orbMulP Flat matrix of orbital contributions [nspin x ntotal_orbitals]
 * @param state   Constraint state (indexing maps)
 * @return Nested vector [nspin][iat][iw]
 */
std::vector<std::vector<std::vector<double>>> convert_orbital_matrix(
    const ModuleBase::matrix& orbMulP,
    const ScState& state);

/**
 * @brief Calculate magnetic moments from converted orbital matrix.
 *
 * @par Algorithm (nspin=2):
 *   atom_mag = sum(orbMulP[0][iat]) - sum(orbMulP[1][iat]); Mi[iat].z = atom_mag
 *
 * @par Algorithm (nspin=4):
 *   total_charge_soc[1..3] = Tr(rho * sigma_{x,y,z}) -> Mi x/y/z.
 *   Components below sc_thr_ are set to 0.0 to avoid noise.
 *
 * @param AorbMulP Nested vector [nspin][iat][iw] from convert_orbital_matrix()
 * @param state    Constraint state (Mi_ written)
 */
void calculate_mw_from_orbitals(const std::vector<std::vector<std::vector<double>>>& AorbMulP,
                                ScState& state);

/**
 * @brief Accumulate magnetic moment contributions from mu*density matrix.
 *
 * @details For distributed matrices (ScaLAPACK), only the local processor's
 * elements are accumulated. The ParaV mapping converts global indices to
 * local row/column indices.
 *
 * @par nspin=4 spinor decomposition
 * The mud matrix stores the 2x2 spinor blocks interleaved:
 *   Global index 2j -> spin-up component
 *   Global index 2j+1 -> spin-down component
 * The Pauli matrix traces are:
 *   M0 (charge): mud(k1,k1).real + mud(k2,k2).real
 *   M3 (Mz):     mud(k1,k1).real - mud(k2,k2).real
 *   M1 (Mx):     mud(k1,k2).real + mud(k2,k1).real
 *   M2 (My):    -mud(k1,k2).imag + mud(k2,k1).imag
 *
 * @param MecMulP Output matrix [4 x nw/2]: MecMulP[0]=charge, [1]=Mx, [2]=My, [3]=Mz
 * @param mud     Input mu*density matrix (column-major)
 * @param nw      Total number of orbitals
 * @param isk     Spin index (0 or 1 for nspin=2)
 * @param state   Constraint state (nspin_)
 * @param pv      Parallel orbitals distribution mapping
 */
void collect_mw(ModuleBase::matrix& MecMulP,
                const ModuleBase::ComplexMatrix& mud,
                int nw,
                int isk,
                const ScState& state,
                const Parallel_Orbitals* pv);

} // namespace lcao
} // namespace spinconstrain

#endif // __LCAO
#endif // DELTASPIN_LCAO_MI_H
