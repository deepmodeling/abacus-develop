// ============================================================
// DFPT_PW::init + DFPT_PW::Impl::build_occ_kq implementation
// with helper extraction, moved from dfpt_pw.cpp so the driver
// translation unit stays below the coding-rule 500-line budget.
// ============================================================

#include "dfpt_pw_impl.h"

#include "dfpt_kq_basis.h"
#include "dfpt_nscf.h"
#include "dfpt_pert.h"
#include "dfpt_phon.h"
#include "dfpt_pw_data.h"
#include "dfpt_q0.h"
#include "dfpt_rho.h"
#include "source_base/constants.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include <cmath>
#include <cstdlib>
#include <map>
#include <memory>
#include <sstream>
#include <vector>

namespace ModuleDFPT
{

void DFPT_PW::Impl::check_metallic_occ(const ModuleBase::matrix &wg) const
{
    // Metallic-sampling guard: the Sternheimer/projector flow treats every
    // band as either fully occupied or empty and carries no d(mu)/dtau
    // response, so a sampling whose smearing Fermi level cuts a band (wg
    // strictly between 0 and the full reference) yields force constants
    // wrong at the 100% level while still converging cleanly. Reject it
    // explicitly (C4 defers metallic DFPT); negligible gauss tails
    // (relative weight < 1e-3) are tolerated as the insulator limit.
    const double frac_weight_tol = 1.0e-3; ///< empirical parameter: relative band weight treated as metallic
    for (int ik = 0; ik < wg.nr; ++ik)
    {
        const double wref = wg(ik, 0);
        if (wref <= 0.0)
        {
            continue;
        }
        for (int ib = 0; ib < wg.nc; ++ib)
        {
            const double rel = wg(ik, ib) / wref;
            if (rel > frac_weight_tol && rel < 1.0 - frac_weight_tol)
            {
                std::stringstream msg;
                msg << "fractional band occupation at (ik=" << ik << ", ib=" << ib << ", wg=" << wg(ik, ib)
                    << "): metallic DFPT (smearing occupations crossing the"
                       " Fermi level) is not supported; reduce smearing sigma"
                       " or use an insulating k sampling.";
                ModuleBase::WARNING_QUIT("DFPT_PW::init", msg.str());
            }
        }
    }
}

void DFPT_PW::Impl::setup_q_list(UnitCell &ucell)
{
    if (!qfile_.empty())
    {
        qlist_.read_from_file(qfile_, ucell);
        if (qlist_.get_nq() == 0)
        {
            ModuleBase::WARNING_QUIT("DFPT_PW::init", "failed to read the DFPT q-point file: " + qfile_);
        }
    }
    else
    {
        std::vector<int> mp_grid = {nqx_, nqy_, nqz_};
        qlist_.generate_mesh(ucell, ucell.symm, mp_grid, true);
    }
    // The uniform E-field response (q0) pipeline runs at q index 0 and is
    // only defined at Gamma, and the dyn-file polar tensors plus LO-TO
    // are keyed to the Gamma entry the same way. Reject, at input time,
    // the orderings that would silently mis-attach those (Gamma not
    // first, or a duplicate Gamma later); arbitrary orderings stay legal
    // for plain phonon runs (compute_q0 off). See validate_q_ordering.
    validate_q_ordering(qlist_, data_.get_compute_q0());
}

void DFPT_PW::Impl::init_submodules(
    const DFPT_PW_InitContext &ctx,
    int nk,
    int nbands,
    int npw_max,
    int nrxx,
    int nspin,
    int nat)
{
    // plain-mixing coefficient: the response Jacobian has strongly
    // negative eigenvalues concentrated on the smallest-G shells (the
    // Coulomb stiffness 4pi/G^2; measured lambda ~ -2.2 on {111}/{200}
    // for the diamond smoke case), so the coefficient must stay below
    // 2 / (1 + |lambda_min|); the INPUT default 0.4 keeps margin up to
    // |lambda| ~ 3; the alternative is dfpt_mix_type = "rho_kerker", the
    // screen f_g = |G+q|^2 / (|G+q|^2 + a^2) in 1/lat0^2 units (a^2 via
    // dfpt_kerker_a2), which stabilizes those shells at beta up to 1;
    // both arrive through the INPUT setters (dfpt_mix_beta /
    // dfpt_mix_type / dfpt_kerker_a2)
    const double mix_beta = mix_beta_;
    const std::string mix_type = mix_type_;
    const double kerker_a2 = kerker_a2_;
    // DFPT_Rho handles density mixing only (the E-field loop in
    // dfpt_pw_q0.cpp always mixes the density, while the v_sc modes do
    // their own potential mixing in solve_displacement), so pass the
    // density-space counterpart: v_sc_plain -> rho_plain, v_sc_broyden
    // -> rho_broyden (the E-field solve must not silently fall back to
    // the plain beta mix when the INPUT selects a Broyden scheme)
    std::string rho_mix_type = mix_type;
    if (mix_type == "v_sc_plain")
    {
        rho_mix_type = "rho_plain";
    }
    else if (mix_type == "v_sc_broyden")
    {
        rho_mix_type = "rho_broyden";
    }
    rho_.init({nspin, nrxx, ctx.pw_rho, ctx.pw_wfc, ctx.ucell->G, rho_mix_type, mix_beta, kerker_a2, broyden_ndim_});
    phon_.init(*ctx.ucell, ctx.pw_rho, &pert_);
    // nlcc direct-term providers (QE dynmat_nlcc + dynmatcc equivalents):
    // the pure core-response builder and the ground-state XC potential
    phon_.set_rho(&rho_);
    // the Impl-owned copy: the context pointer refers to caller stack
    // state that dies when the esolver init returns
    // always pass the pointer (never nullptr): a rank with an empty
    // local field must still enter accum_nlcc_diag's collective
    // participation; the size checks inside gate only the local work
    phon_.set_vxc_r(&vxc_r_);
    q0_.init(*ctx.ucell, ctx.pw_rho, ctx.pw_wfc, &pert_);
    hamilt_.reset(new DFPT_HamiltShift(*ctx.ucell, ctx.pw_rho, ctx.pw_wfc, *ctx.veff_r, &pert_));
    data_.init(&qlist_, nk, nbands, npw_max, nrxx, nspin, nat, ctx.dftu);
}

void DFPT_PW::init(const DFPT_PW_InitContext &ctx)
{
    ModuleBase::TITLE("DFPT_PW", "init");
    ModuleBase::timer::start("DFPT_PW", "init");
    // ctx.{psi,wg,veff_r,eig} are valid non-null pointers even in skeleton
    // mode: the empty Psi / matrix objects still carry the queryable nk /
    // nbands / nrxx shape fields consumed by the init helpers; ucell is
    // always a non-null per the public wrapper (it takes a reference).
    UnitCell *const ucell = ctx.ucell;
    const psi::Psi<std::complex<double>> *const psi = ctx.psi;
    ModulePW::PW_Basis *const pw_rho = ctx.pw_rho;
    ModulePW::PW_Basis_K *const pw_wfc = ctx.pw_wfc;
    Structure_Factor *const sf = ctx.sf;

    pimpl_->ucell_ = ucell;
    pimpl_->gs_psi_ = *psi;
    pimpl_->pw_rho_ = pw_rho;
    pimpl_->pw_wfc_ = pw_wfc;
    pimpl_->sf_ = sf;
    pimpl_->veff_r_ = *ctx.veff_r;
    pimpl_->vxc_r_.clear();
    if (ctx.vxc_r != nullptr)
    {
        pimpl_->vxc_r_ = *ctx.vxc_r;
    }
    pimpl_->wg_ = *ctx.wg;
    pimpl_->eig_ = *ctx.eig;
    pimpl_->xc_ = ctx.xc;
    pimpl_->nelec_ = ctx.nelec;
    pimpl_->ecutwfc_ = ctx.ecutwfc;
    pimpl_->dftu_ = ctx.dftu;
    pimpl_->myrank_ = ctx.myrank;

    pimpl_->check_metallic_occ(*ctx.wg);

    // DFT+U guard: the ground state now supports PW-basis DFT+U and wires a
    // provider when dft_plus_u is enabled, but every DFPT U hook
    // (DFPT_Rho::cal_docc, DFPT_Pert::build_dv_u, DFPT_Q0 born/docc
    // contractions, DFPT_Phon::dftu_onsite) is a no-op reservation (U0).
    // Running anyway would converge cleanly while silently dropping the
    // whole first-order U response, so reject explicitly until U1 lands
    // (same fail-loud pattern as the metallic-sampling guard above).
    if (ctx.dftu != nullptr)
    {
        ModuleBase::WARNING_QUIT(
            "DFPT_PW::init", "DFT+U with DFPT is not supported yet: the "
                             "first-order U response is not implemented "
                             "(U0 reservation); rerun with dft_plus_u 0.");
    }

    pimpl_->setup_q_list(*ucell);

    const int nk = psi->get_nk();
    const int nbands = psi->get_nbands();
    const int npw_max = psi->get_current_ngk();
    const int nrxx = (pw_rho != nullptr) ? pw_rho->nrxx : 0;
    const int nspin = 1;
    const int nat = ucell->nat;

    if (pw_rho != nullptr && pw_wfc != nullptr && sf != nullptr)
    {
        pimpl_->pert_.init(*ucell, pw_rho, pw_wfc, *sf);
        pimpl_->init_submodules(ctx, nk, nbands, npw_max, nrxx, nspin, nat);
    }
    else
    {
        pimpl_->phon_.init(*ucell, nullptr, nullptr);
        pimpl_->data_.init(&pimpl_->qlist_, nk, nbands, npw_max, nrxx, nspin, nat, ctx.dftu);
    }
    ModuleBase::timer::end("DFPT_PW", "init");
}

int DFPT_PW::Impl::occupied_band_count(int ik) const
{
    // insulator: the occupied-band count is k-independent; the metallic
    // guard in init() has already rejected a cut Fermi surface
    int nocc = 0;
    const int nbands = wg_.nc;

    for (int m = 0; m < nbands; ++m)
    {
        if (dfpt_band_occupied(wg_, ik, m))
        {

            nocc = m + 1;
        }
    }
    return nocc;
}

void DFPT_PW::Impl::build_occ_kq(int q_idx)
{
    ModuleBase::TITLE("DFPT_PW", "build_occ_kq");
    ModuleBase::timer::start("DFPT_PW", "build_occ_kq");
    const int nk = pw_wfc_->nks;
    occ_kq_.assign(nk, std::vector<std::vector<std::complex<double>>>());
    eig_kq_.assign(nk, std::vector<double>());
    const ModuleBase::Vector3<double> q_frac = data_.get_qvec(q_idx);
    const ModuleBase::Vector3<double> q_cart = q_frac * ucell_->G;
    const int nocc = occupied_band_count(0);

    // QE alignment: the k+q occupied states come from an nscf
    // diagonalization of the full H(k+q) on the k+q ball (run_nscf ->
    // non_scf_ph), not from any copy of the ground-state list. The
    // block-Davidson runs on DFPT_HamiltShift with the shift zeroed; the
    // per-band Sternheimer shift keeps using the ground-state eigenvalue
    // eig_(ik, ib) of the coupled k band, so eig_ stays the source there.
    const int lin_max = 200;          ///< empirical parameter: Davidson outer iterations per k+q ball
    const double eig_tol = 1.0e-5;    ///< empirical parameter: residual norm target (Ry); occ_kq_ feeds the
                                      ///< P_c projector only, whose leakage at this level stays far below
                                      ///< the response convergence threshold
    const double precond_floor = 0.1; ///< empirical parameter: kinetic-preconditioner floor (1/lat0^2); small
                                      ///< enough that the near-Gamma plane waves (which dominate the
                                      ///< valence top) are not suppressed in the expansion directions
    const double tpiba2 = ucell_->tpiba2;
    // Gamma q: H(k+q) = H_GS(k) identically, so the k+q occupied states
    // ARE the ground-state states -- the diagonalization would only add
    // its own eigensolver noise (and a cross-toolchain variance at that).
    // QE treats the Gamma q the same way (run_nscf's lgamma_iq branch
    // skips the nscf re-diagonalization). Take the occupied states from
    // the ground-state wavefunctions directly; q != 0 keeps the nscf path.
    const bool q_is_gamma = (q_cart.norm2() < 1.0e-12);
    const ModuleBase::Matrix3 ginv = pw_wfc_->G.Inverse();
    for (int ik = 0; ik < nk; ++ik)
    {
        if (q_is_gamma)
        {
            hamilt_->set_context(q_cart, ik);
            hamilt_->set_shift(0.0);
            const DFPT_KQ_Basis &kq = hamilt_->basis();
            const int npwk = kq.get_npwk();
            // ground-state per-k index by the reciprocal-basis integer
            // triple (pw_wfc_ stores a per-k ball layout: read it through
            // getgcar)
            std::map<std::vector<int>, int> jgl_of_n;
            for (int jgl = 0; jgl < pw_wfc_->npwk[ik]; ++jgl)
            {
                const ModuleBase::Vector3<double> gf = pw_wfc_->getgcar(ik, jgl) * ginv;
                const std::vector<int> key = {static_cast<int>(std::round(gf.x)), static_cast<int>(std::round(gf.y)),
                    static_cast<int>(std::round(gf.z))};
                jgl_of_n[key] = jgl;
            }
            occ_kq_[ik].clear();
            eig_kq_[ik].clear();
            for (int m = 0; m < gs_psi_.get_nbands(); ++m)
            {
                if (!dfpt_band_occupied(wg_, ik, m))
                {
                    continue;
                }
                std::vector<std::complex<double>> state(npwk, std::complex<double>(0.0, 0.0));
                for (int igl = 0; igl < npwk; ++igl)
                {
                    const ModuleBase::Vector3<double> gf = kq.get_gcar(igl) * ginv;
                    const std::vector<int> key = {static_cast<int>(std::round(gf.x)),
                        static_cast<int>(std::round(gf.y)), static_cast<int>(std::round(gf.z))};
                    const std::map<std::vector<int>, int>::const_iterator it = jgl_of_n.find(key);
                    if (it != jgl_of_n.end())
                    {
                        state[igl] = gs_psi_(ik, m, it->second);
                    }
                }
                occ_kq_[ik].push_back(std::move(state));
                eig_kq_[ik].push_back(eig_(ik, m));
            }
            continue;
        }
        hamilt_->set_context(q_cart, ik);
        hamilt_->set_shift(0.0);
        const DFPT_KQ_Basis &kq = hamilt_->basis();
        const int npwk = kq.get_npwk();
        std::vector<double> precond(npwk);
        for (int igl = 0; igl < npwk; ++igl)
        {
            precond[igl] = kq.get_gk2(igl) * tpiba2 + precond_floor;
        }
        std::vector<std::vector<std::complex<double>>> occ_ik;
        std::vector<double> eig_ik;
        const int it = DFPT_NSCF::lowest_eigpairs(*hamilt_, precond, nocc, lin_max, eig_tol, occ_ik, eig_ik);
        if (it != 0)
        {
            std::ostringstream oss;
            oss << "the nscf diagonalization at k index " << ik << " did not reach the eigenvalue "
                << "residual target within " << lin_max << " outer iterations (last pass " << it
                << "); raise dfpt_conv_thr precision handling or inspect the spectrum";
            ModuleBase::WARNING_QUIT("DFPT_PW::build_occ_kq", oss.str());
        }
        occ_kq_[ik] = std::move(occ_ik);
        eig_kq_[ik] = std::move(eig_ik);
    }
    last_q_ = q_idx;
    last_ik_ = -1;
    ModuleBase::timer::end("DFPT_PW", "build_occ_kq");
}

} // namespace ModuleDFPT
