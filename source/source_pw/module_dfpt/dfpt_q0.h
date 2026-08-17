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
 *            * Re[<u_v|r_a|u_c><u_c|r_b|u_v>] / (eps_c - eps_v) / Nk
 * Born charges from dP/dtau (King-Smith/Resta Berry phases; the m sum runs
 * over ALL bands, occupied and empty, m != v):
 *   Z*_k,ab = Z_k delta_ab - (4/Nk) sum_{k,v occ,m!=v} wg
 *             * Re[<u_v|dV/dtau_{k,b}|u_m><u_m|r_a|u_v>] / (eps_m - eps_v)
 * The bare displacement potential dV/dtau comes from DFPT_Pert (C1) at
 * q = 0; the absolute calibration of both expressions is pinned by the
 * diamond end-to-end test in C7 (structure/symmetry by the C6 tests).
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
    
private:
    UnitCell* ucell_ = nullptr;
    ModulePW::PW_Basis* pw_rho_ = nullptr;
    ModulePW::PW_Basis_K* pw_wfc_ = nullptr;
    DFPT_Pert* pert_ = nullptr;
};

} // namespace ModuleDFPT

#endif // DFPT_Q0_H
