#ifndef DFPT_PW_IMPL_H
#define DFPT_PW_IMPL_H

#include "dfpt_hamilt_shift.h"
#include "dfpt_kq_basis.h"
#include "dfpt_metal.h"
#include "dfpt_pert.h"
#include "dfpt_phon.h"
#include "dfpt_pw.h"
#include "dfpt_pw_data.h"
#include "dfpt_q0.h"
#include "dfpt_rho.h"
#include "dfpt_stern.h"
#include "source_base/matrix.h"
#include "source_base/vector3.h"
#include "source_cell/qlist.h"
#include "source_cell/unitcell.h"
#include "source_psi/psi.h"

#include <memory>
#include <string>
#include <vector>

namespace ModulePW
{
class PW_Basis;
class PW_Basis_K;
} // namespace ModulePW

class Plus_U_Base;
class Structure_Factor;

namespace ModuleDFPT
{

class XC_First_Order;

/// Private implementation body of DFPT_PW (opaque in the public header).
///
/// Only the constructors / destructor and the outward-facing per-solve
/// hooks are public; all data members and per-solve helpers are private.
class DFPT_PW::Impl
{
  public:
    Impl();
    ~Impl();

    /// occupied-state projector set at k+q for every k of this q
    /// (commensurate q: kvec_d[ik] + q must be a k point of the
    /// ground-state list mod lattice)
    void build_occ_kq(int q_idx);

    /// one self-consistent Sternheimer cycle for the displacement
    /// (iat, idir) at q; returns the achieved density residual (zero
    /// when unwired)
    double solve_displacement(int q_idx, int iat, int idir);

    /// position legs Y^a_{k,v} = P_c x_a|psi_{k,v}> of the q = 0 mesh
    void solve_pos_resp(int q_idx);

    /// E-field SCF response dpsi^E,a of the q = 0 mesh
    void solve_efield_resp(int q_idx);

    /// per-iteration assembly of the screened response potential
    /// v_sc^q = v_Hartree(drho_in) + xc(drho_in) on the shared grid;
    /// shared by solve_displacement and solve_efield_resp so the
    /// real-space Hartree+XC add loop is not duplicated
    void assemble_v_sc(const ModuleBase::Vector3<double>& q_cart,
                       const std::vector<std::complex<double>>& drho_in_g,
                       std::vector<std::complex<double>>& v_sc_r) const;

  private:
    friend class DFPT_PW;

    // ----- initialisation helpers (dfpt_pw_init.cpp) -----
    void check_metallic_occ(const ModuleBase::matrix& wg) const;
    void setup_q_list(UnitCell& ucell);
    void init_submodules(const DFPT_PW_InitContext& ctx,
                         int nk,
                         int nbands,
                         int npw_max,
                         int nrxx,
                         int nspin,
                         int nat);

    // ----- build_occ_kq helpers (dfpt_pw_init.cpp) -----
    int match_commensurate_kq(int ik,
                              const ModuleBase::Vector3<double>& q_frac,
                              double tol,
                              ModuleBase::Vector3<int>& dn_out) const;
    void copy_occ_state_ball(int ik,
                             int ikq,
                             const ModuleBase::Vector3<int>& dn,
                             const DFPT_KQ_Basis& kq,
                             const ModuleBase::Matrix3& ginv,
                             std::vector<std::vector<std::complex<double>>>& occ_ik) const;

    // ----- run() dispatch helpers (dfpt_pw_run.cpp) -----
    void run_q0_pre(int q_idx);
    double run_displacement_irrep_pass(int q_idx, int irrep);
    void run_q0_post(int q_idx);
    void run_assemble(int q_idx);

    // ----- solve_displacement helpers (dfpt_pw_solve.cpp) -----
    double sternheimer_per_band(int ik,
                                int ib,
                                const std::vector<std::vector<std::complex<double>>>& dv_sc,
                                int nbands,
                                int lin_max,
                                double lin_thr);
    void stash_converged_disp_response(int q_idx,
                                       int iat,
                                       int idir,
                                       const std::vector<std::complex<double>>& v_sc_r_last,
                                       int nk,
                                       int nbands);

    // ----- q=0 solve helpers (dfpt_pw_q0.cpp) -----
    void vel_diag_part(int ik, int a, int nbands, std::vector<std::vector<std::complex<double>>>& vel) const;
    void vel_nl_per_atom(int ik,
                         int a,
                         int it,
                         int ia,
                         int nbands,
                         const std::vector<ModuleBase::Vector3<double>>& gk,
                         std::vector<std::vector<std::complex<double>>>& vel);
    void pos_per_band_solve(int ik,
                            int a,
                            int nbands,
                            int lin_max,
                            double lin_thr,
                            std::vector<std::vector<std::vector<std::complex<double>>>>& yvec);
    void efield_per_band_solve(int ik,
                               int a,
                               int nbands,
                               int lin_max,
                               double lin_thr,
                               const std::vector<std::vector<std::vector<std::complex<double>>>>& yr,
                               const std::vector<std::vector<std::complex<double>>>& dv_sc);
    void stash_dpsi_efield(int q_idx, int a, int nk, int nbands);

    bool wired() const;

    DFPT_PW_Data data_;
    DFPT_Pert pert_;
    DFPT_Stern stern_;
    DFPT_Rho rho_;
    DFPT_Phon phon_;
    DFPT_Q0 q0_;
    DFPT_Metal metal_;
    ModuleCell::QList qlist_;
    std::unique_ptr<DFPT_HamiltShift> hamilt_;

    psi::Psi<std::complex<double>> gs_psi_;
    UnitCell* ucell_ = nullptr;
    ModulePW::PW_Basis* pw_rho_ = nullptr;
    ModulePW::PW_Basis_K* pw_wfc_ = nullptr;
    Structure_Factor* sf_ = nullptr;
    std::vector<double> veff_r_;
    ModuleBase::matrix wg_;
    ModuleBase::matrix eig_;
    const XC_First_Order* xc_ = nullptr;
    double nelec_ = 0.0;
    double ecutwfc_ = 0.0;
    const Plus_U_Base* dftu_ = nullptr;

    ///< occupied states at k+q on the k+q G list, [ik][occ m][igl];
    ///< rebuilt per q (they depend on q and k only)
    std::vector<std::vector<std::vector<std::complex<double>>>> occ_kq_;
    ///< remembers the (q_idx, ik) the shifted operator was last cached at
    int last_q_ = -1;
    int last_ik_ = -1;
    std::vector<int> ikq_of_k_;

    int nqx_ = 1, nqy_ = 1, nqz_ = 1;
    std::string qfile_;
    double conv_thr_ = 1e-8;
    int max_iter_ = 100;
    double mix_beta_ = 0.4;
};

} // namespace ModuleDFPT

#endif // DFPT_PW_IMPL_H
