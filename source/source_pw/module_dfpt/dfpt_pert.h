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

    /// C5: read access to the ground-state wfc basis for the dynamical-matrix
    /// contractions in DFPT_Phon::accumulate_electron.
    ModulePW::PW_Basis_K* get_pw_wfc() const { return pw_wfc_; }
    
    void build_dv(int q_idx, int atom_idx, int dir, DFPT_PW_Data& data);
    
    void apply_dv(int q_idx, int k_idx, const psi::Psi<std::complex<double>>& psi, 
                  DFPT_PW_Data& data);
    
    void build_efield(const ModuleBase::Vector3<double>& field, DFPT_PW_Data& data);

    /// C5: real-space kernel of the same-atom second-order LOCAL potential
    /// d2V_loc(r) = d^2 V_loc / d tau_{da} d tau_{db} (both displacements on
    /// the SAME atom; each derivative contributes i w_dir, w = Delta + q, so
    /// the reciprocal kernel is -w_da w_db Vloc(|w|) exp(i w.tau)). Returned
    /// on the shared real-space grid; its expectation value with |u(r)|^2
    /// enters the electronic dynamical matrix (anharmonic term).
    void d2vloc_r(int atom_idx, int da, int db,
                  const ModuleBase::Vector3<double>& q_cart,
                  std::vector<std::complex<double>>& dv2_r) const;

    /// C5: same-atom second-order NONLOCAL potential acting on psi,
    /// chi_n(G'') = (d^2 Vnl / d tau_{da} d tau_{db}) |psi_n> on the k+q
    /// basis (normal-conserving separable case). The four terms come from the
    /// phase derivatives of the out (k+q+G'') and in (k+G') projectors of the
    /// SAME displaced atom; they reduce to zero for a uniform translation at
    /// q=0 (acoustic consistency).
    void apply_d2vnl(int atom_idx, int da, int db,
                     const ModuleBase::Vector3<double>& q_cart,
                     const psi::Psi<std::complex<double>>& psi, int k_idx,
                     std::vector<std::vector<std::complex<double>>>& d2v_psi) const;

    /// Build the beta-projector array (in the ABACUS vkb convention) for a
    /// single atom on an arbitrary k-shifted reciprocal vector list:
    ///   vkb[mu][ielem] = (-i)^l * Ylm(Ghat) * (4pi/sqrt(Omega) *
    ///                        integral beta(r) j_l(g r) r dr) * exp(i G.tau)
    /// with G in 2*pi/lat0 units, g = |G| * tpiba (bohr^-1), tau in bohr.
    /// Usable for both the incoming k basis (G = k+G') and the outgoing DFPT
    /// k+q basis (G = k+q+G''), so the atomic phase is correct on either side.
    /// Public since C6: DFPT_Q0 reuses it for the velocity operator.
    void build_vkb(int it, int ia,
                   const std::vector<ModuleBase::Vector3<double>>& gk,
                   std::vector<std::vector<std::complex<double>>>& vkb) const;

    /// C6: analytic derivative of the beta projector with respect to the
    /// list shift k_dir (the same shift build_vkb is evaluated at), three
    /// terms: the atomic phase (i 2pi tau_dir), the radial chain rule
    /// (vq'(g) * tpiba * Ghat_dir, central finite difference of radial_vq)
    /// and the real-harmonic direction derivative (grad_real_ylm chain
    /// (e_dir - ghat ghat_dir)/|G|). Feeds the dV_nl/dk part of the
    /// velocity operator in DFPT_Q0::pos_matrix.
    void build_vkb_dk(int it, int ia, int dir,
                      const std::vector<ModuleBase::Vector3<double>>& gk,
                      std::vector<std::vector<std::complex<double>>>& vkb,
                      std::vector<std::vector<std::complex<double>>>& dvkb) const;

    /// radial part (4pi/sqrt(Omega)) Integral beta(r) j_l(g r) r dr at g (bohr^-1)
    double radial_vq(int it, int ib, double g) const;

    /// C7: apply a complex real-space potential on the shared FFT grid to
    /// every band of psi (k basis), delivering |v psi> on the k+q basis.
    /// The potential is the q-shifted complex periodic amplitude (the same
    /// convention as dv_rc); the DFPT self-consistent loop uses it for the
    /// screened response potential (Hartree + XC) of the mixed density.
    void apply_vr(int q_idx, int k_idx,
                  const std::vector<std::complex<double>>& v_rc,
                  const psi::Psi<std::complex<double>>& psi,
                  const ModuleBase::Vector3<double>& q_cart,
                  std::vector<std::vector<std::complex<double>>>& dv_psi) const;

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
    ///                       * exp(-i (Delta+q).tau_atom) * ...
    /// (GS structure-factor convention exp(-2pi g.tau); the sign/coefficient
    /// is the exact derivative of the local potential with respect to the
    /// atomic displacement).
    void dVloc_dtau(int atom_idx, int dir, const ModuleBase::Vector3<double>& q, 
                    std::vector<std::complex<double>>& dv);
    
    /// C1: first-order NONLOCAL potential acting on psi (normal-conserving
    /// separable case), for one displaced atom in direction dir.
    /// Uses the identity (GS exp(-2pi gk.tau) projector convention)
    ///   dVnl/dtau_a |psi> = -i (k+q+G'')_a * (Vnl |psi>)
    ///                       + Vnl[ i (k+G')_a |psi> ]
    /// so only two applications of the ground-state nonlocal operator on the
    /// DFPT k+q outgoing basis are needed (dsVnl contribution per pair is
    /// i (q+G''-G')_a times the zero-order matrix element).
    /// USPP/ultrasoft and spin-orbit projectors are rejected for now.
    void dVnl_dtau(int atom_idx, int dir, const ModuleBase::Vector3<double>& q,
                   const psi::Psi<std::complex<double>>& psi, int k_idx,
                    std::vector<std::vector<std::complex<double>>>& dv_psi);

    /// real spherical harmonic Y_{l,m}(g_hat), orthonormal convention, l<=2.
    double real_ylm(int l, int m, const ModuleBase::Vector3<double>& ghat) const;
    /// gradient of real_ylm with respect to the unit vector ghat, l<=2
    /// (dY/dghat returned per cartesian component).
    void grad_real_ylm(int l, int m, const ModuleBase::Vector3<double>& ghat,
                       double grad[3]) const;

    /// General (nonlocal and local) part of apply_dv for the compartments that
    /// live in real space (local potential); the |psi> product requires the
    /// shared real-space grid of pw_rho_/pw_wfc_.
    void real_space_dv(int q_idx, int k_idx,
                       const psi::Psi<std::complex<double>>& psi,
                       DFPT_PW_Data& data,
                       const DFPT_KQ_Basis& kq,
                       std::vector<std::vector<std::complex<double>>>& dv_psi) const;

    /// shared core of real_space_dv / apply_vr: phase-free cyclic convolution
    /// of v_rc with every band of psi, scattered/gathered between the k+q
    /// list and the rho grid through the FFT-cell triple.
    void apply_vr_core(int k_idx,
                       const std::vector<std::complex<double>>& v_rc,
                       const psi::Psi<std::complex<double>>& psi,
                       const DFPT_KQ_Basis& kq,
                       std::vector<std::vector<std::complex<double>>>& dv_psi) const;

    /// first-order Hubbard potential dV_U (U0 reservation, C1 frozen term).
    void build_dv_u(int q_idx, int atom_idx, int dir, DFPT_PW_Data& data);
};

} // namespace ModuleDFPT

#endif // DFPT_PERT_H