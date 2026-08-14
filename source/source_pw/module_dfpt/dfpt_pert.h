// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in the design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================

#ifndef DFPT_PERT_H
#define DFPT_PERT_H

#include "dfpt_kq_basis.h"
#include "dfpt_pw_data.h"
#include "source_cell/unitcell.h"
#include "source_psi/psi.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_basis/module_pw/pw_basis_k.h"

class Structure_Factor;

namespace ModuleDFPT {

class DFPT_Pert {
public:
    DFPT_Pert();
    ~DFPT_Pert();
    
    void init(UnitCell& ucell, ModulePW::PW_Basis* pw_rho, 
              ModulePW::PW_Basis_K* pw_wfc, Structure_Factor& sf);
    
    void build_dv(int q_idx, int atom_idx, int dir, DFPT_PW_Data& data);
    
    void apply_dv(int q_idx, int k_idx, const psi::Psi<std::complex<double>>& psi, 
                  DFPT_PW_Data& data);
    
    void build_efield(const ModuleBase::Vector3<double>& field, DFPT_PW_Data& data);

private:
    UnitCell* ucell_ = nullptr;
    ModulePW::PW_Basis* pw_rho_ = nullptr;
    ModulePW::PW_Basis_K* pw_wfc_ = nullptr;
    Structure_Factor* sf_ = nullptr;

    /// C1: first-order LOCAL potential dVloc_dtau (per displaced atom).
    /// Grid helper: reconstruct the cartesian reciprocal vector (in 2*pi/lat0
    /// units) of rho-grid index ig from the shared FFT-grid (ix,iy,iz) mapping.
    void rho_gvec(int ig, ModuleBase::Vector3<double>& gcar) const;
    /// The local pseudopotential Vloc(g^2) at an arbitrary magnitude:
    ///  Coulomb atoms use the analytic form, numeric pseudopotentials reuse the
    ///  radial-mesh Fourier transform of vl_pw.cpp::vloc_of_g at |g| themselves.
    double vloc_at_g(int it, double g2) const;
    /// linear atom index -> (type, picture) of ucell_.
    void atom_index(int atom_idx, int& it, int& ia) const;

    /// First-order asymmetric-part local potential on the rho grid:
    ///   dVloc_dtau(Delta) = -i (Delta+q).direction * Vloc(|Delta+q|)
    ///                       * exp(i (Delta+q).tau_atom) * ...
    /// The sign/coefficient convention is the exact derivative of the local
    /// potential with respect to the atomic displacement (see source; the unit
    /// test checks it against a finite difference of the full potential incl.
    /// the atomic phase).
    void dVloc_dtau(int atom_idx, int dir, const ModuleBase::Vector3<double>& q, 
                    std::vector<std::complex<double>>& dv);
    
    /// C1: first-order NONLOCAL potential acting on psi (normal-conserving
    /// separable case), for one displaced atom in direction dir.
    /// Uses the identity
    ///   dVnl/dtau_a |psi> = i (k+q+G'')_a * (Vnl |psi>) - Vnl[ i (k+G')_a |psi> ]
    /// so only two applications of the ground-state nonlocal operator on the
    /// DFPT k+q outgoing basis are needed (dsVnl contribution per pair is
    /// i (q+G''-G')_a times the zero-order matrix element).
    /// USPP/ultrasoft and spin-orbit projectors are rejected for now.
    void dVnl_dtau(int atom_idx, int dir, const ModuleBase::Vector3<double>& q,
                   const psi::Psi<std::complex<double>>& psi, int k_idx,
                   std::vector<std::vector<std::complex<double>>>& dv_psi);

    /// Build the beta-projector array (in the ABACUS vkb convention) for a
    /// single atom on an arbitrary k-shifted reciprocal vector list:
    ///   vkb[mu][ielem] = (-i)^l * Ylm(Ghat) * (4pi/sqrt(Omega) *
    ///                        integral beta(r) j_l(g r) r dr) * exp(i G.tau)
    /// with G in 2*pi/lat0 units, g = |G| * tpiba (bohr^-1), tau in bohr.
    /// Usable for both the incoming k basis (G = k+G') and the outgoing DFPT
    /// k+q basis (G = k+q+G''), so the atomic phase is correct on either side.
    void build_vkb(int it, int ia,
                   const std::vector<ModuleBase::Vector3<double>>& gk,
                   std::vector<std::vector<std::complex<double>>>& vkb) const;
    /// radial part (4pi/sqrt(Omega)) Integral beta(r) j_l(g r) r dr at g (bohr^-1)
    double radial_vq(int it, int ib, double g) const;
    /// real spherical harmonic Y_{l,m}(g_hat), orthonormal convention, l<=2.
    double real_ylm(int l, int m, const ModuleBase::Vector3<double>& ghat) const;

    /// General (nonlocal and local) part of apply_dv for the compartments that
    /// live in real space (local potential); the |psi> product requires the
    /// shared real-space grid of pw_rho_/pw_wfc_.
    void real_space_dv(int q_idx, int k_idx,
                       const psi::Psi<std::complex<double>>& psi,
                       DFPT_PW_Data& data,
                       const DFPT_KQ_Basis& kq,
                       std::vector<std::vector<std::complex<double>>>& dv_psi) const;

    /// first-order Hubbard potential dV_U (U0 reservation, C1 frozen term).
    void build_dv_u(int q_idx, int atom_idx, int dir, DFPT_PW_Data& data);
};

} // namespace ModuleDFPT

#endif // DFPT_PERT_H