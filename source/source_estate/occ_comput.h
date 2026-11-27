#ifndef OCC_COMPUT_H
#define OCC_COMPUT_H

#include <complex>

namespace elecstate
{

/**
 * @brief Accumulate per-projector 2x2 occupation blocks from the onsite
 *        projector coefficients for a single k-point.
 *
 * For each projector (l,m) (global index iprj), the occupation block is
 *
 *   rho^{ss'}_{iprj} = sum_i w_{k,i} * conj(proj^s_{i,iprj}) * proj^{s'}_{i,iprj}
 *
 * where proj^s_{i,iprj} = <alpha_{iprj}|Psi^s_{k,i}> and w_{k,i} is the band
 * occupation weight. The block is stored as
 *   occ_block[iprj*4 + {0,1,2,3}] = {rho^{up,up}, rho^{up,dn},
 *                                    rho^{dn,up}, rho^{dn,dn}}.
 *
 * The 2x2 layout is a storage convention, not an assumption of spin physics;
 * the content of the four slots depends on nspin:
 *   nspin=1: rho is split evenly into occ[0]/occ[3], off-diagonal slots stay
 *            zero. No Pauli structure is involved; the even split only makes
 *            the magnetization readout M ~ occ[0] - occ[3] vanish.
 *   nspin=2: isk=0 accumulates into occ[0], isk=1 into occ[3]; the two
 *            diagonal slots are simply two independent spin channels.
 *   nspin=4: both spinor components cross, filling all four slots with a
 *            genuine 2x2 Hermitian spin density matrix (the only case where
 *            the Pauli decomposition rho = (n*I + M.sigma)/2 applies).
 *
 * Accumulation semantics (+=): the caller zeroes occ_block before the
 * k-point loop and performs the MPI reduction after it.
 *
 * @param proj      onsite projector coefficients <alpha|psi>, laid out as
 *                  (nbands*npol) x nkb with spinor components offset by nkb
 * @param wg_ik     band occupation weights for this k-point (nbands values)
 * @param nbands    number of bands
 * @param npol      number of spinor components (1 collinear, 2 non-collinear)
 * @param nkb       total number of projectors
 * @param nspin     1, 2 or 4
 * @param isk       spin channel of this k-point (0/1); used only when nspin=2
 * @param nh_iat    number of projectors per atom (nat values)
 * @param nat       number of atoms
 * @param occ_block [in,out] blocks to accumulate into, size sum(nh_iat)*4
 */
void occ_from_proj(
    const std::complex<double>* proj,
    const double* wg_ik,
    const int nbands,
    const int npol,
    const int nkb,
    const int nspin,
    const int isk,
    const int* nh_iat,
    const int nat,
    std::complex<double>* occ_block);

} // namespace elecstate

#endif // OCC_COMPUT_H
