// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in the design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================

#ifndef DFPT_Q0_H
#define DFPT_Q0_H

#include "dfpt_pw_data.h"
#include "source_cell/unitcell.h"
#include "source_psi/psi.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_basis/module_pw/pw_basis_k.h"

namespace ModuleDFPT {

class DFPT_Pert;

/**
 * @brief q -> 0 response: dielectric tensor, Born charges, LO-TO (C6).
 *
 * The position operator is ill-defined for periodic states, so the
 * periodic-gauge matrix elements are obtained through the velocity
 * (commutator) form (Gonze & Lee, PRB 55, 10355 (1997)), m != n:
 *   <u_m|r_dir|u_n> = -i <u_m|dH/dk_dir|u_n> / (tpiba (eps_m - eps_n)),
 *   dH/dk_dir = 2 tpiba^2 (k+G)_dir (diagonal kinetic part)
 *             + dV_nl/dk_dir (build_vkb_dk; V_loc is k-independent),
 * with the k derivative in dimensionless 2*pi/lat0 units (matching
 * build_vkb_dk) so r comes out in bohr; [r, V_U] of DFT+U is a
 * documented U0 reservation (the onsite projector is nonlocal, so the
 * commutator does not vanish with U on).
 *
 * Dielectric tensor (insulating, ABACUS Ry a.u. with wg carrying the spin
 * degeneracy; the extra 1/(eps_c - eps_v) is the length-gauge denominator,
 * consistent with the oscillator-strength sum rule):
 *   eps_ab = delta_ab + (8 pi / Omega) sum_{k,v occ,c emp} wg
 *            * Re[<u_v|r_a|u_c><u_c|r_b|u_v>] / (eps_c - eps_v)
 * Born charges from dP/dtau, the screened displacement leg paired with the
 * SOLVED conduction-projected position response (QE zstar_eu/add_zstar_ue
 * anchoring; Gonze-Lee screened form). The position leg
 *   Y^a_{k,v} = P_c x_a|psi_{k,v}>,  (H(k)-eps_v) Y = P_c [H,x_a]|psi_v>,
   * with the commutator rhs [H,x_a]|psi> = -(i/tpiba) dH/dk_a|psi> (the
   * same velocity operator as above), is solved exactly by Sternheimer
   * solves in DFPT_PW (solve_pos_resp, stashed per direction in the shared
   * data) and therefore carries the complete conduction-space response; the
   * eigenvector-truncated r-matrix contraction of the du form is only its
   * nbands-cut approximation. With dpsi^kappa(scf) the converged q = 0
   * Sternheimer displacement responses:
   *   Z*_k,ab = Z_k delta_ab - 2 sum_{k,v occ} wg
   *             * Re <dpsi^{k,b}_scf,v | Y^a_v>
   * (wg carries the spin degeneracy, so the prefactor is the -2*wk of
   * add_zstar_ue). By the symmetry of the mixed second derivative of the
   * total energy this equals the transposed leg
   * -2*sum wg*Re<dpsi^E_a(scf)|dV_ext^{k,b}|psi_v> that QE's zstar_eu
   * computes with the electric-field responses; only one leg is needed.
   * The dpsi^kappa Sternheimer gauge (<psi_occ|dpsi> = 0) drops the
   * occupied-occupied block of x exactly. The diamond C7 target
   * (Z* -> 0 by inversion + ASR) requires the screened dpsi.
 * With a symmetry-reduced k list both sums run over the irreducible k and
 * each partial tensor chi(k) is star-averaged: the physical partial at a
 * rotated star member Rk is R chi(k) R^T, and atom-resolved (Born) partials
 * are credited to the image atom under R. With symmetry off the stored list
 * is the full mesh and the star machinery degenerates to the identity.
 * The absolute calibration of both expressions is pinned by the diamond
 * end-to-end test in C7 (structure/symmetry by the C6 tests).
 */
class DFPT_Q0 {
public:
    DFPT_Q0();
    ~DFPT_Q0();

    void init(UnitCell& ucell, ModulePW::PW_Basis* pw_rho,
              ModulePW::PW_Basis_K* pw_wfc, DFPT_Pert* pert);

    void compute_eps(const psi::Psi<std::complex<double>>& psi,
                     const ModuleBase::matrix& wg,
                     const ModuleBase::matrix& eig, DFPT_PW_Data& data);

    void compute_born(const psi::Psi<std::complex<double>>& psi,
                      const ModuleBase::matrix& wg,
                      const ModuleBase::matrix& eig, DFPT_PW_Data& data);

    void compute_q0_response(DFPT_PW_Data& data);

    /// position-operator matrix elements r_mat[ik][m][n].d =
    /// <u_{m,k}|r_d|u_{n,k}> (m != n), periodic gauge (velocity form);
    /// eig is the ground-state eigenvalue matrix (nk x nbands, Ry).
    void pos_matrix(const psi::Psi<std::complex<double>>& psi,
                    const ModuleBase::matrix& eig,
                    std::vector<std::vector<std::vector<ModuleBase::Vector3<std::complex<double>>>>>& r_mat);

    // ---- k-star rotation of the symmetry-reduced q=0 tensor sums ----
    // One entry per DISTINCT folded star member Rk of an irreducible k:
    // the representative operation in cartesian COLUMN form (rotate_tensor
    // applies chi' = R chi R^T directly) plus the atom map iat -> image
    // atom under the same operation (built from the direct space
    // gmatrix/gtrans pair; species map to themselves).
    struct StarMember {
        ModuleBase::Matrix3 cart;      ///< defaults to the identity
        std::vector<int> atom_map;     ///< empty means the identity map
    };
    std::vector<std::vector<StarMember>> stars_;  ///< [ik] -> star members

    /// rebuild stars_ for the stored k list (nk points); falls back to a
    /// single identity member per k when the point group is unavailable
    /// (symmetry off / unreduced mesh) or an atom map fails
    void build_stars(int nk);

    /// chi_rot(a,b) = sum_{a'b'} R(a,a') R(b,b') chi(a',b') of a 3x3
    /// partial tensor under a cartesian rotation
    static void rotate_tensor(const ModuleBase::Matrix3& r,
                              const ModuleBase::matrix& chi,
                              double (&chi_rot)[9]);

private:
    UnitCell* ucell_ = nullptr;
    ModulePW::PW_Basis* pw_rho_ = nullptr;
    ModulePW::PW_Basis_K* pw_wfc_ = nullptr;
    DFPT_Pert* pert_ = nullptr;
};

} // namespace ModuleDFPT

#endif // DFPT_Q0_H
