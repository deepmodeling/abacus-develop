// ============================================================
// DFPT_PW::init + DFPT_PW::Impl::build_occ_kq implementation
// with helper extraction, moved from dfpt_pw.cpp so the driver
// translation unit stays below the coding-rule 500-line budget.
// ============================================================

#include "dfpt_pw_impl.h"

#include "dfpt_kq_basis.h"
#include "dfpt_pert.h"
#include "dfpt_phon.h"
#include "dfpt_pw_data.h"
#include "dfpt_q0.h"
#include "dfpt_rho.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"

#include <cstdlib>
#include <map>
#include <memory>
#include <sstream>
#include <vector>

namespace ModuleDFPT
{

void DFPT_PW::Impl::check_metallic_occ(const ModuleBase::matrix& wg) const
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

void DFPT_PW::Impl::setup_q_list(UnitCell& ucell)
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
}

void DFPT_PW::Impl::init_submodules(const DFPT_PW_InitContext& ctx,
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
    // |lambda| ~ 3; the alternative is mix_type = "kerker", the screen
    // f_g = |G+q|^2 / (|G+q|^2 + a^2) in 1/lat0^2 units (a^2 via
    // DFPT_KERKER_A2), which stabilizes those shells at beta up to 1;
    // the env knobs are design-phase calibration aids
    double mix_beta = mix_beta_;
    if (const char* env_beta = getenv("DFPT_MIX_BETA"))
    {
        const double parsed = atof(env_beta);
        if (parsed > 0.0 && parsed <= 1.0)
        {
            mix_beta = parsed;
        }
    }
    std::string mix_type = "plain";
    if (const char* env_type = getenv("DFPT_MIX_TYPE"))
    {
        const std::string parsed = env_type;
        if (parsed == "plain" || parsed == "kerker")
        {
            mix_type = parsed;
        }
    }
    double kerker_a2 = 1.0;
    if (const char* env_a2 = getenv("DFPT_KERKER_A2"))
    {
        const double parsed = atof(env_a2);
        if (parsed > 0.0)
        {
            kerker_a2 = parsed;
        }
    }
    rho_.init({nspin, nrxx, ctx.pw_rho, ctx.pw_wfc, ctx.ucell->G, mix_type, mix_beta, kerker_a2});
    phon_.init(*ctx.ucell, ctx.pw_rho, &pert_);
    q0_.init(*ctx.ucell, ctx.pw_rho, ctx.pw_wfc, &pert_);
    hamilt_.reset(new DFPT_HamiltShift(*ctx.ucell, ctx.pw_rho, ctx.pw_wfc, *ctx.veff_r, &pert_));
    data_.init(&qlist_, nk, nbands, npw_max, nrxx, nspin, nat, ctx.dftu);
}

void DFPT_PW::init(const DFPT_PW_InitContext& ctx)
{
    ModuleBase::TITLE("DFPT_PW", "init");
    ModuleBase::timer::start("DFPT_PW", "init");
    // ctx.{psi,wg,veff_r,eig} are valid non-null pointers even in skeleton
    // mode: the empty Psi / matrix objects still carry the queryable nk /
    // nbands / nrxx shape fields consumed by the init helpers; ucell is
    // always a non-null per the public wrapper (it takes a reference).
    UnitCell* const ucell = ctx.ucell;
    const psi::Psi<std::complex<double>>* const psi = ctx.psi;
    ModulePW::PW_Basis* const pw_rho = ctx.pw_rho;
    ModulePW::PW_Basis_K* const pw_wfc = ctx.pw_wfc;
    Structure_Factor* const sf = ctx.sf;

    pimpl_->ucell_ = ucell;
    pimpl_->gs_psi_ = *psi;
    pimpl_->pw_rho_ = pw_rho;
    pimpl_->pw_wfc_ = pw_wfc;
    pimpl_->sf_ = sf;
    pimpl_->veff_r_ = *ctx.veff_r;
    pimpl_->wg_ = *ctx.wg;
    pimpl_->eig_ = *ctx.eig;
    pimpl_->xc_ = ctx.xc;
    pimpl_->nelec_ = ctx.nelec;
    pimpl_->ecutwfc_ = ctx.ecutwfc;
    pimpl_->dftu_ = ctx.dftu;

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
        ModuleBase::WARNING_QUIT("DFPT_PW::init",
                                 "DFT+U with DFPT is not supported yet: the "
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

int DFPT_PW::Impl::match_commensurate_kq(int ik,
                                         const ModuleBase::Vector3<double>& q_frac,
                                         double tol,
                                         ModuleBase::Vector3<int>& dn_out) const
{
    // k+q folded into [0,1) direct coordinates must be a ground-state k
    // point (DFPT q meshes are commensurate with the k mesh)
    const ModuleBase::Vector3<double> target = pw_wfc_->kvec_d[ik] + q_frac;
    const int nk = pw_wfc_->nks;
    for (int j = 0; j < nk; ++j)
    {
        const ModuleBase::Vector3<double>& kj = pw_wfc_->kvec_d[j];
        const double rx = std::round(kj.x - target.x);
        const double ry = std::round(kj.y - target.y);
        const double rz = std::round(kj.z - target.z);
        if (std::abs(kj.x - target.x - rx) < tol && std::abs(kj.y - target.y - ry) < tol
            && std::abs(kj.z - target.z - rz) < tol)
        {
            dn_out.x = static_cast<int>(rx);
            dn_out.y = static_cast<int>(ry);
            dn_out.z = static_cast<int>(rz);
            return j;
        }
    }
    std::ostringstream oss;
    oss << "k+q is not a point of the ground-state k list: the DFPT "
           "q mesh must be commensurate with the k mesh (and inside "
           "the first Brillouin zone). ik="
        << ik << " k_d=(" << pw_wfc_->kvec_d[ik].x << "," << pw_wfc_->kvec_d[ik].y << ","
        << pw_wfc_->kvec_d[ik].z << ") q_d=(" << q_frac.x << "," << q_frac.y << "," << q_frac.z << ") k+q=("
        << target.x << "," << target.y << "," << target.z << ") nk=" << nk;
    ModuleBase::WARNING_QUIT("DFPT_PW::build_occ_kq", oss.str());
}

void DFPT_PW::Impl::copy_occ_state_ball(int ik,
                                        int ikq,
                                        const ModuleBase::Vector3<int>& dn,
                                        const DFPT_KQ_Basis& kq,
                                        const ModuleBase::Matrix3& ginv,
                                        std::vector<std::vector<std::complex<double>>>& occ_ik) const
{
    // reciprocal-basis integer triple -> per-k index of the ikq ball
    // (pw_wfc_ is a PW_Basis_K whose gcar holds a per-k ball layout,
    // not the parent-class global-ig layout: read it through getgcar)
    std::map<std::vector<int>, int> jgl_of_n;
    for (int jgl = 0; jgl < pw_wfc_->npwk[ikq]; ++jgl)
    {
        const ModuleBase::Vector3<double> gf = pw_wfc_->getgcar(ikq, jgl) * ginv;
        const std::vector<int> key = {static_cast<int>(std::round(gf.x)),
                                      static_cast<int>(std::round(gf.y)),
                                      static_cast<int>(std::round(gf.z))};
        jgl_of_n[key] = jgl;
    }

    const int npw_kq = kq.get_npwk();
    const int nbands = gs_psi_.get_nbands();
    for (int m = 0; m < nbands; ++m)
    {
        if (!dfpt_band_occupied(wg_, ikq, m))
        {
            continue; // empty at k+q: outside the P_c projector
        }
        std::vector<std::complex<double>> state(npw_kq, std::complex<double>(0.0, 0.0));
        for (int igl = 0; igl < npw_kq; ++igl)
        {
            const ModuleBase::Vector3<double> gf = kq.get_gcar(igl) * ginv;
            const std::vector<int> key = {static_cast<int>(std::round(gf.x)) + dn.x,
                                          static_cast<int>(std::round(gf.y)) + dn.y,
                                          static_cast<int>(std::round(gf.z)) + dn.z};
            const auto it = jgl_of_n.find(key);
            if (it != jgl_of_n.end())
            {
                state[igl] = gs_psi_(ikq, m, it->second);
            }
        }
        occ_ik.push_back(std::move(state));
    }
}

void DFPT_PW::Impl::build_occ_kq(int q_idx)
{
    ModuleBase::TITLE("DFPT_PW", "build_occ_kq");
    ModuleBase::timer::start("DFPT_PW", "build_occ_kq");
    const int nk = pw_wfc_->nks;
    occ_kq_.assign(nk, std::vector<std::vector<std::complex<double>>>());
    ikq_of_k_.assign(nk, -1);
    const ModuleBase::Vector3<double> q_frac = data_.get_qvec(q_idx);
    const ModuleBase::Vector3<double> q_cart = q_frac * ucell_->G;
    const double kmatch_tol = 1.0e-6; ///< empirical parameter: folded fractional k-list match tolerance
    for (int ik = 0; ik < nk; ++ik)
    {
        ModuleBase::Vector3<int> dn(0, 0, 0);
        const int ikq = match_commensurate_kq(ik, q_frac, kmatch_tol, dn);
        ikq_of_k_[ik] = ikq;

        DFPT_KQ_Basis kq;
        kq.init(pw_wfc_, pw_rho_, q_cart, ik);
        const ModuleBase::Matrix3 ginv = pw_wfc_->G.Inverse();
        copy_occ_state_ball(ik, ikq, dn, kq, ginv, occ_kq_[ik]);
    }
    last_q_ = q_idx;
    last_ik_ = -1;
    ModuleBase::timer::end("DFPT_PW", "build_occ_kq");
}

} // namespace ModuleDFPT
