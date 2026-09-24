#ifndef DM_FROM_PSI_H
#define DM_FROM_PSI_H

#include "density_matrix.h"
#include "source_base/matrix.h"
#include "source_psi/psi.h"

namespace module_dm
{
/**
 * @brief Calculate the k-point density matrix from wavefunctions for the Gamma-only case
 *
 * For every k-point block of the wavefunction this computes
 *     dmk(iw1, iw2) = sum_ib wg(ik, ib) * wfc(ib, iw1) * wfc(ib, iw2),
 * i.e. the wavefunction weights are real so no conjugation is needed.
 *
 * @param ParaV 2D block-cyclic distribution descriptor of the orbitals
 * @param wg band weights, wg(ik, ib_global) for the global band index
 * @param wfc wavefunction coefficients, shape (nk, nbands_local, nbasis_local)
 * @param DM output density matrix, one block per k-point
 */
void dm_from_psi(const Parallel_Orbitals* ParaV,
                const ModuleBase::matrix& wg,
                const psi::Psi<double>& wfc,
                DensityMatrix<double, double>& DM);

/**
 * @brief Calculate the k-point density matrix from wavefunctions for the multi-k case
 *
 * For every k-point block of the wavefunction this computes
 *     dmk(iw1, iw2) = sum_ib wg(ik, ib) * conj(wfc(ib, iw1)) * wfc(ib, iw2).
 * The conjugation is applied to the FIRST basis index ("conj-first" storage):
 * the stored block is the transpose of the physical one-body density matrix
 * P = C diag(wg) C^H. cal_dmr(), SOC magnetization and Mulliken analysis all
 * rely on this convention, see dm_tools.cpp and
 * unittests/test_soc_magnetization_roundtrip.cpp.
 *
 * @tparam TR real type of the real-space density matrix of DM
 * @param ParaV 2D block-cyclic distribution descriptor of the orbitals
 * @param wg band weights, wg(ik, ib_global) for the global band index
 * @param wfc wavefunction coefficients, shape (nk, nbands_local, nbasis_local)
 * @param DM output density matrix, one block per k-point
 */
template <typename TR>
void dm_from_psi(const Parallel_Orbitals* ParaV,
                const ModuleBase::matrix& wg,
                const psi::Psi<std::complex<double>>& wfc,
                DensityMatrix<std::complex<double>, TR>& DM);

/**
 * @brief Calculate one k-point block of the density matrix (Gamma-only case)
 *
 * Single-k-point worker behind dm_from_psi(); exposed for callers that own the
 * output storage themselves. The conjugation/weighting contract is the same as
 * dm_from_psi(), the caller only provides the output slot.
 *
 * @param ParaV 2D block-cyclic distribution descriptor of the orbitals
 * @param wg band weights, wg(ik, ib_global) for the global band index
 * @param ik k-point index of the block to calculate
 * @param wfc wavefunction coefficients
 * @param dmk_out output block, distributed according to ParaV->desc (serial
 *                builds: a dense nbasis_local x nbasis_local column-major matrix)
 */
void dmk_from_psi(const Parallel_Orbitals* ParaV,
                 const ModuleBase::matrix& wg,
                 const int ik,
                 const psi::Psi<double>& wfc,
                 double* dmk_out);

/**
 * @brief Calculate one k-point block of the density matrix (multi-k case)
 *
 * Single-k-point worker behind dm_from_psi(), using the conj-first convention
 * dmk(iw1, iw2) = sum_ib wg(ik, ib) * conj(wfc(ib, iw1)) * wfc(ib, iw2).
 *
 * @param ParaV 2D block-cyclic distribution descriptor of the orbitals
 * @param wg band weights, wg(ik, ib_global) for the global band index
 * @param ik k-point index of the block to calculate
 * @param wfc wavefunction coefficients
 * @param dmk_out output block, distributed according to ParaV->desc (serial
 *                builds: a dense nbasis_local x nbasis_local column-major matrix)
 */
void dmk_from_psi(const Parallel_Orbitals* ParaV,
                 const ModuleBase::matrix& wg,
                 const int ik,
                 const psi::Psi<std::complex<double>>& wfc,
                 std::complex<double>* dmk_out);
} // namespace module_dm
#endif
