#ifndef MI_TOOLS_H
#define MI_TOOLS_H

#include <complex>
#include <vector>

#include "source_base/vector3.h"

/**
 * @file mi_tools.h
 * @brief Free-function utilities for computing atomic magnetic moments (Mi)
 *        from intermediate quantities (e.g. becp projector coefficients)
 *        in the DeltaSpin module.
 *
 * @par Rationale
 * These functions do not depend on the SpinConstrain<TK> template state and
 * are factored out as free functions to:
 *   - Reduce coupling to the singleton class members
 *   - Allow reuse by both LCAO and PW basis paths
 *   - Be independently unit-testable
 */

namespace spinconstrain
{

/**
 * @brief Convert spinor occupation matrix to magnetic moment vector using Pauli matrices.
 *
 * @details For a two-component spinor wavefunction, the spin density matrix is:
 *   rho = |a|^2    a*b  |   = | (1+Mz)/2    (Mx-iMy)/2 |
 *         |b*a    |b|^2  |     | (Mx+iMy)/2   (1-Mz)/2  |
 * The magnetic moment components are extracted via Pauli matrix traces:
 *   Mx = Tr(rho * sigma_x) = occ[1] + occ[2]           (real part)
 *   My = Tr(rho * sigma_y) = -Im(occ[1] - occ[2])      (from sigma_y = [[0,-i],[i,0]])
 *   Mz = Tr(rho * sigma_z) = occ[0] - occ[3]            (real part)
 * where occ = {|a|^2, a*b, b*a, |b|^2} from becp coefficients.
 *
 * @param occ 4-element array of occupation matrix elements (complex)
 * @param weight k-point weight for integration
 * @return 3D magnetic moment vector (Mx, My, Mz) in Bohr magnetons
 */
inline ModuleBase::Vector3<double> pauli_to_moment(const std::complex<double> occ[4], double weight)
{
    return ModuleBase::Vector3<double>(
        weight * (occ[1] + occ[2]).real(),
        weight * (occ[1] - occ[2]).imag(),
        weight * (occ[0] - occ[3]).real()
    );
}

/**
 * @brief Accumulate atomic magnetic moments from becp coefficients for one k-point.
 *
 * @details For npol=2 (nspin=4), computes full Pauli decomposition:
 *   occ[0] = sum(becp_up^* * becp_up), occ[1] = sum(becp_up^* * becp_dn),
 *   occ[2] = sum(becp_dn^* * becp_up), occ[3] = sum(becp_dn^* * becp_dn)
 *   Mi = pauli_to_moment(occ, weight)
 * For npol=1 (nspin=2), only z-component:
 *   occ = sum(|becp|^2), Mi.z += weight * occ * spin_sign
 *
 * @param becp Projector coefficients <alpha_{l,m}|psi_{k,i}>
 * @param nkb Total number of projectors
 * @param nbands Number of bands
 * @param npol Number of spinor components (1 for collinear, 2 for non-collinear)
 * @param spin_sign +1 for spin-up, -1 for spin-down (nspin=2 only); unused for npol=2
 * @param wg_ik Band occupation weights for this k-point (from Fermi-Dirac)
 * @param nh_iat Array of projector counts per atom: nh_iat[iat] = nproj for atom iat
 * @param mi [in,out] Magnetic moments vector to accumulate into (size = nat)
 */
void accumulate_Mi_from_becp(const std::complex<double>* becp,
                            int nkb,
                            int nbands,
                            int npol,
                            int spin_sign,
                            const double* wg_ik,
                            const int* nh_iat,
                            std::vector<ModuleBase::Vector3<double>>& mi);

} // namespace spinconstrain

#endif // MI_TOOLS_H
