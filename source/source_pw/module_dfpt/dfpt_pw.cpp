// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================

#include "dfpt_pw.h"
#include "dfpt_pw_data.h"
#include "dfpt_pert.h"
#include "dfpt_stern.h"
#include "dfpt_rho.h"
#include "dfpt_phon.h"
#include "dfpt_q0.h"
#include "dfpt_metal.h"
#include "dfpt_hamilt_shift.h"
#include "dfpt_kq_basis.h"
#include "source_base/constants.h"
#include "source_base/global_function.h"
#include <sstream>
#include "source_cell/qlist.h"
#include "source_pw/module_pwdft/stru_fac.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <complex>
#include <vector>

namespace ModuleDFPT {

class DFPT_PW::Impl {
public:
    Impl() {}
    ~Impl()
    {
        delete hamilt_;
    }

    DFPT_PW_Data data_;
    DFPT_Pert pert_;
    DFPT_Stern stern_;
    DFPT_Rho rho_;
    DFPT_Phon phon_;
    DFPT_Q0 q0_;
    DFPT_Metal metal_;
    ModuleCell::QList qlist_;
    DFPT_HamiltShift* hamilt_ = nullptr;

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
    const Plus_U* dftu_ = nullptr;

    ///< occupied states at k+q on the k+q G list, [ik][occ m][igl];
    /// rebuilt per q (they depend on q and k only)
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

    bool wired() const { return pw_rho_ != nullptr && pw_wfc_ != nullptr; }

    /// occupied-state projector set at k+q for every k of this q (commensurate
    /// q: kvec_d[ik] + q must be a k point of the ground-state list mod lattice)
    void build_occ_kq(int q_idx);

    /// one self-consistent Sternheimer cycle for the displacement (iat, idir)
    /// at q; returns the achieved density residual (zero when unwired)
    double solve_displacement(int q_idx, int iat, int idir);

    /// position legs Y^a_{k,v} = P_c x_a|psi_{k,v}> of the q = 0 mesh
    /// (velocity-rhs Sternheimer solves, one per direction; stashed through
    /// data as the exact position leg of the screened Born charges)
    void solve_pos_resp(int q_idx);

    /// E-field SCF response dpsi^E,a of the q = 0 mesh (QE solve_e +
    /// dfpt_kernel form: fixed point on the rhs -(Y^a + dV_sc^E,a|psi>)
    /// with the screening assembly of solve_displacement)
    void solve_efield_resp(int q_idx);
};

DFPT_PW::DFPT_PW() : pimpl_(new Impl()) {}

DFPT_PW::~DFPT_PW() {
    delete pimpl_;
}

void DFPT_PW::init(UnitCell& ucell, const psi::Psi<std::complex<double>>& psi,
                   ModulePW::PW_Basis* pw_rho, ModulePW::PW_Basis_K* pw_wfc,
                   Structure_Factor* sf, const std::vector<double>& veff_r,
                   const ModuleBase::matrix& wg, const ModuleBase::matrix& eig,
                   const XC_First_Order* xc,
                   double nelec, double ecutwfc, const Plus_U* dftu) {
    pimpl_->ucell_ = &ucell;
    pimpl_->gs_psi_ = psi;
    pimpl_->pw_rho_ = pw_rho;
    pimpl_->pw_wfc_ = pw_wfc;
    pimpl_->sf_ = sf;
    pimpl_->veff_r_ = veff_r;
    pimpl_->wg_ = wg;
    pimpl_->eig_ = eig;

    // Metallic-sampling guard: the Sternheimer/projector flow treats every
    // band as either fully occupied or empty and carries no d(mu)/dtau
    // response, so a sampling whose smearing Fermi level cuts a band (wg
    // strictly between 0 and the full reference) yields force constants
    // wrong at the 100% level while still converging cleanly. Reject it
    // explicitly (C4 defers metallic DFPT); negligible gauss tails
    // (relative weight < 1e-3) are tolerated as the insulator limit.
    for (int ik = 0; ik < wg.nr; ++ik) {
        const double wref = wg(ik, 0);
        if (wref <= 0.0) {
            continue;
        }
        for (int ib = 0; ib < wg.nc; ++ib) {
            const double rel = wg(ik, ib) / wref;
            if (rel > 1.0e-3 && rel < 1.0 - 1.0e-3) {
                std::stringstream msg;
                msg << "fractional band occupation at (ik=" << ik
                    << ", ib=" << ib << ", wg=" << wg(ik, ib)
                    << "): metallic DFPT (smearing occupations crossing the"
                       " Fermi level) is not supported; reduce smearing sigma"
                       " or use an insulating k sampling.";
                ModuleBase::WARNING_QUIT("DFPT_PW::init", msg.str());
            }
        }
    }
    pimpl_->xc_ = xc;
    pimpl_->nelec_ = nelec;
    pimpl_->ecutwfc_ = ecutwfc;
    pimpl_->dftu_ = dftu;

    // q points: an explicit q list file overrides the Monkhorst-Pack mesh
    if (!pimpl_->qfile_.empty()) {
        pimpl_->qlist_.read_from_file(pimpl_->qfile_, ucell);
        if (pimpl_->qlist_.get_nq() == 0) {
            ModuleBase::WARNING_QUIT("DFPT_PW::init",
                                     "failed to read the DFPT q-point file: " + pimpl_->qfile_);
        }
    } else {
        std::vector<int> mp_grid = {pimpl_->nqx_, pimpl_->nqy_, pimpl_->nqz_};
        pimpl_->qlist_.generate_mesh(ucell, ucell.symm, mp_grid, true);
    }

    int nq = pimpl_->qlist_.get_nq();
    int nk = psi.get_nk();
    int nbands = psi.get_nbands();
    int npw_max = psi.get_current_ngk();
    int nrxx = (pw_rho != nullptr) ? pw_rho->nrxx : 0;
    int nspin = 1;
    int nat = ucell.nat;

    if (pw_rho != nullptr && pw_wfc != nullptr && sf != nullptr) {
        pimpl_->pert_.init(ucell, pw_rho, pw_wfc, *sf);
        // plain-mixing coefficient: the response Jacobian has strongly
        // negative eigenvalues concentrated on the smallest-G shells (the
        // Coulomb stiffness 4pi/G^2; measured lambda ~ -2.2 on {111}/{200}
        // for the diamond smoke case), so the coefficient must stay below
        // 2 / (1 + |lambda_min|); the INPUT default 0.4 keeps margin up to
        // |lambda| ~ 3; the alternative is mix_type = "kerker", the screen
        // f_g = |G+q|^2 / (|G+q|^2 + a^2) in 1/lat0^2 units (a^2 via
        // DFPT_KERKER_A2), which stabilizes those shells at beta up to 1;
        // the env knobs are design-phase calibration aids
        double mix_beta = pimpl_->mix_beta_;
        if (const char* env_beta = getenv("DFPT_MIX_BETA")) {
            const double parsed = atof(env_beta);
            if (parsed > 0.0 && parsed <= 1.0) {
                mix_beta = parsed;
            }
        }
        std::string mix_type = "plain";
        if (const char* env_type = getenv("DFPT_MIX_TYPE")) {
            const std::string parsed = env_type;
            if (parsed == "plain" || parsed == "kerker") {
                mix_type = parsed;
            }
        }
        double kerker_a2 = 1.0;
        if (const char* env_a2 = getenv("DFPT_KERKER_A2")) {
            const double parsed = atof(env_a2);
            if (parsed > 0.0) {
                kerker_a2 = parsed;
            }
        }
        pimpl_->rho_.init(nspin, nrxx, pw_rho, pw_wfc, ucell.G, mix_type, mix_beta, kerker_a2);
        pimpl_->phon_.init(ucell, pw_rho, &pimpl_->pert_);
        pimpl_->q0_.init(ucell, pw_rho, pw_wfc, &pimpl_->pert_);
        delete pimpl_->hamilt_;
        pimpl_->hamilt_ = new DFPT_HamiltShift(ucell, pw_rho, pw_wfc, veff_r, &pimpl_->pert_);
    } else {
        pimpl_->phon_.init(ucell, nullptr, nullptr);
    }
    pimpl_->data_.init(&pimpl_->qlist_, nk, nbands, npw_max, nrxx, nspin, nat, dftu);
}

bool DFPT_PW::get_with_u() const {
    return pimpl_->data_.with_u();
}

bool DFPT_PW::get_u_active() const {
    return pimpl_->data_.u_active();
}

void DFPT_PW::Impl::build_occ_kq(int q_idx) {
    const int nk = pw_wfc_->nks;
    occ_kq_.assign(nk, std::vector<std::vector<std::complex<double>>>());
    ikq_of_k_.assign(nk, -1);
    const ModuleBase::Vector3<double> q_frac = data_.get_qvec(q_idx);
    const ModuleBase::Vector3<double> q_cart = q_frac * ucell_->G;
    for (int ik = 0; ik < nk; ++ik) {
        // k+q folded into [0,1) direct coordinates must be a ground-state k
        // point (DFPT q meshes are commensurate with the k mesh)
        const ModuleBase::Vector3<double> target = pw_wfc_->kvec_d[ik] + q_frac;
        int ikq = -1;
        for (int j = 0; j < nk; ++j) {
            const ModuleBase::Vector3<double>& kj = pw_wfc_->kvec_d[j];
            const double rx = std::round(kj.x - target.x);
            const double ry = std::round(kj.y - target.y);
            const double rz = std::round(kj.z - target.z);
            if (std::abs(kj.x - target.x - rx) < 1.0e-6
                && std::abs(kj.y - target.y - ry) < 1.0e-6
                && std::abs(kj.z - target.z - rz) < 1.0e-6) {
                ikq = j;
                break;
            }
        }
        if (ikq < 0) {
            std::ostringstream oss;
            oss << "k+q is not a point of the ground-state k list: the DFPT "
                   "q mesh must be commensurate with the k mesh (and inside "
                   "the first Brillouin zone). ik=" << ik
                << " k_d=(" << pw_wfc_->kvec_d[ik].x << "," << pw_wfc_->kvec_d[ik].y
                << "," << pw_wfc_->kvec_d[ik].z << ") q_d=(" << q_frac.x << ","
                << q_frac.y << "," << q_frac.z << ") k+q=(" << target.x << ","
                << target.y << "," << target.z << ") nk=" << nk;
            ModuleBase::WARNING_QUIT("DFPT_PW::build_occ_kq", oss.str());
        }
        ikq_of_k_[ik] = ikq;

        DFPT_KQ_Basis kq;
        kq.init(pw_wfc_, pw_rho_, q_cart, ik);
        const int npw_kq = kq.get_npwk();

        // The congruence match above may fold k+q onto a *different label*
        // of the same physical point (e.g. k lists holding both (1/2,0,0)
        // and (-1/2,0,0), which differ by a reciprocal lattice vector).
        // The two balls then enumerate different G labels: a state of the
        // ikq ball with vector G' coincides physically with the k+q-ball
        // vector G when G' + k(ijq) == G + k(ik) + q, i.e.
        // G' = G + dn with dn = k_d(ik) + q - k_d(ikq) integer in
        // reciprocal-basis coordinates. Coincident FFT cells identify the
        // same G only for dn = 0, so match through the G vectors instead.
        const ModuleBase::Vector3<double> dn = pw_wfc_->kvec_d[ik] + q_frac
                                               - pw_wfc_->kvec_d[ikq];
        const double dnr[3] = {std::round(dn.x), std::round(dn.y), std::round(dn.z)};
        if (std::abs(dn.x - dnr[0]) > 1.0e-6 || std::abs(dn.y - dnr[1]) > 1.0e-6
            || std::abs(dn.z - dnr[2]) > 1.0e-6) {
            ModuleBase::WARNING_QUIT("DFPT_PW::build_occ_kq",
                                     "k+q folds onto a k-list entry with a "
                                     "non-integer reciprocal offset.");
        }
        const int dn_i[3] = {static_cast<int>(dnr[0]),
                             static_cast<int>(dnr[1]),
                             static_cast<int>(dnr[2])};
        const ModuleBase::Matrix3 ginv = pw_wfc_->G.Inverse();
        // reciprocal-basis integer triple -> per-k index of the ikq ball
        // (pw_wfc_ is a PW_Basis_K whose gcar holds a per-k ball layout,
        // not the parent-class global-ig layout: read it through getgcar)
        std::map<std::vector<int>, int> jgl_of_n;
        for (int jgl = 0; jgl < pw_wfc_->npwk[ikq]; ++jgl) {
            const ModuleBase::Vector3<double> gf
                = pw_wfc_->getgcar(ikq, jgl) * ginv;
            const std::vector<int> key = {static_cast<int>(std::round(gf.x)),
                                          static_cast<int>(std::round(gf.y)),
                                          static_cast<int>(std::round(gf.z))};
            jgl_of_n[key] = jgl;
        }

        const int nbands = gs_psi_.get_nbands();
        for (int m = 0; m < nbands; ++m) {
            if (!dfpt_band_occupied(wg_, ikq, m)) {
                continue; // empty at k+q: outside the P_c projector
            }
            std::vector<std::complex<double>> state(npw_kq, std::complex<double>(0.0, 0.0));
            for (int igl = 0; igl < npw_kq; ++igl) {
                const ModuleBase::Vector3<double> gf = kq.get_gcar(igl) * ginv;
                const std::vector<int> key
                    = {static_cast<int>(std::round(gf.x)) + dn_i[0],
                       static_cast<int>(std::round(gf.y)) + dn_i[1],
                       static_cast<int>(std::round(gf.z)) + dn_i[2]};
                const auto it = jgl_of_n.find(key);
                if (it != jgl_of_n.end()) {
                    state[igl] = gs_psi_(ikq, m, it->second);
                }
            }
            occ_kq_[ik].push_back(std::move(state));
        }
    }
    last_q_ = q_idx;
    last_ik_ = -1;
}

double DFPT_PW::Impl::solve_displacement(int q_idx, int iat, int idir) {
    if (!wired() || hamilt_ == nullptr) {
        return 0.0;
    }
    const ModuleBase::Vector3<double> q_frac = data_.get_qvec(q_idx);
    const ModuleBase::Vector3<double> q_cart = q_frac * ucell_->G;
    const int nrxx = pw_rho_->nrxx;
    const int nk = gs_psi_.get_nk();
    const int nbands = gs_psi_.get_nbands();

    pert_.build_dv(q_idx, iat, idir, data_);
    rho_.reset_mixing(q_idx);
    // the previous perturbation's stored response must not leak into the
    // first iteration of this one
    data_.set_drho_g(q_idx, 0,
                     std::vector<std::complex<double>>(pw_rho_->npw,
                                                       std::complex<double>(0.0, 0.0)));

    const int lin_max = data_.get_max_iter();
    const double lin_thr = data_.get_conv_thr();

    bool converged = false;
    double residual = 0.0;
    const bool dbg = (getenv("DFPT_DEBUG") != nullptr);
    // last screened response potential (hoisted out of the loop: the 2n+1
    // accumulation below needs the converged v_sc of this displacement)
    std::vector<std::complex<double>> v_sc_r_last;
    for (int iter = 0; iter < max_iter_ && !converged; ++iter) {
        // the per-displacement SCF state (iter / residual / converged) is
        // local to this solve: the DFPT_PW_Data ledger is the per-(q,irrep)
        // outer-pass record kept by run(), and the final residual is
        // returned to the caller for that aggregation (B4)

        // ---- 1. screened response potential from the mixed input density:
        // q-shifted complex periodic amplitude on the shared grid, i.e. the
        // same convention as dv_rc (v_hartree_q acts on the q-shifted
        // coefficients; the XC kernel responds to Re/Im of the amplitude)
        std::vector<std::complex<double>> v_sc_r(nrxx, std::complex<double>(0.0, 0.0));
        const std::vector<std::complex<double>> drho_in_g = data_.get_drho_g(q_idx, 0);
        if (!drho_in_g.empty() && static_cast<int>(drho_in_g.size()) == pw_rho_->npw) {
            std::vector<std::complex<double>> dv_ha_g;
            rho_.v_hartree_q(q_cart, drho_in_g, dv_ha_g);
            std::vector<std::complex<double>> vh_r(nrxx);
            pw_rho_->recip2real(dv_ha_g.data(), vh_r.data());
            for (int ir = 0; ir < nrxx; ++ir) {
                v_sc_r[ir] = vh_r[ir];
            }
            if (xc_ != nullptr) {
                std::vector<std::complex<double>> a_r(nrxx);
                pw_rho_->recip2real(drho_in_g.data(), a_r.data());
                std::vector<std::complex<double>> b_r;
                xc_->apply(a_r, b_r);
                if (static_cast<int>(b_r.size()) == nrxx) {
                    for (int ir = 0; ir < nrxx; ++ir) {
                        v_sc_r[ir] += b_r[ir];
                    }
                }
            }
            if (dbg) {
                double dh = 0.0;
                double dv = 0.0;
                for (int ig = 0; ig < pw_rho_->npw; ++ig) {
                    dh += std::norm(drho_in_g[ig]);
                }
                for (int ir = 0; ir < nrxx; ++ir) {
                    dv += std::norm(v_sc_r[ir]);
                }
                std::cout << "DBG iter=" << iter << " |drho_in_g|=" << std::sqrt(dh)
                          << " |v_sc_r|=" << std::sqrt(dv) << std::endl;
            }
        }
        v_sc_r_last = v_sc_r;

        // ---- 2. Sternheimer solve of every occupied (k, band)
        for (int ik = 0; ik < nk; ++ik) {
            if (static_cast<int>(occ_kq_.size()) <= ik || occ_kq_[ik].empty()) {
                if (dbg) { std::cout << "DBG skip ik=" << ik << " no occ_kq" << std::endl; }
                continue; // no occupied states at k+q: nothing to solve
            }
            // dV_ext |psi_n> for all bands (dVloc convolution + dVnl_dtau)
            pert_.apply_dv(q_idx, ik, gs_psi_, data_);
            // screened response part |v_sc psi_n>
            std::vector<std::vector<std::complex<double>>> dv_sc;
            pert_.apply_vr(q_idx, ik, v_sc_r, gs_psi_, q_cart, dv_sc);
            if (ik != last_ik_ || last_q_ != q_idx) {
                hamilt_->set_context(q_cart, ik);
                last_ik_ = ik;
                if (dbg) {
                    std::cout << "DBG occ_kq nstates=" << occ_kq_[ik].size() << std::endl;
                    for (size_t m = 0; m < occ_kq_[ik].size(); ++m) {
                        double nrm = 0.0;
                        for (size_t i = 0; i < occ_kq_[ik][m].size(); ++i) {
                            nrm += std::norm(occ_kq_[ik][m][i]);
                        }
                        std::cout << "DBG   occ[" << m << "] |psi|^2=" << nrm << std::endl;
                    }
                    // kernel consistency: <psi_m|H(k+q)|psi_m> must equal
                    // eig(ikq, m); the eigenvalue used by set_shift below is
                    // the k-side one (equal only when H is assembled right)
                    for (size_t m = 0; m < occ_kq_[ik].size(); ++m) {
                        hamilt_->set_shift(0.0);
                        std::vector<std::complex<double>> hp(occ_kq_[ik][m].size());
                        hamilt_->apply(occ_kq_[ik][m].data(), hp.data());
                        std::complex<double> dot(0.0, 0.0);
                        for (size_t i = 0; i < hp.size(); ++i) {
                            dot += std::conj(occ_kq_[ik][m][i]) * hp[i];
                        }
                        std::cout << "DBG   <psi_" << m << "|H|psi_" << m << "> = "
                                  << dot.real() << " + i " << dot.imag()
                                  << "  (GS eig " << eig_(ikq_of_k_[ik], static_cast<int>(m)) << ")" << std::endl;
                        std::cout << "DBG   <psi_" << m << "|T+Vnl|psi_" << m << "> = "
                                  << hamilt_->debug_t_vnl(occ_kq_[ik][m]) << std::endl;
                        std::cout << "DBG   <psi_" << m << "|V(wfc-path)|psi_" << m << "> = "
                                  << hamilt_->debug_v_wfc(occ_kq_[ik][m]) << std::endl;
                    }
                }
            }
            for (int ib = 0; ib < nbands; ++ib) {
                if (!dfpt_band_occupied(wg_, ik, ib)) {
                    continue; // unoccupied: no Sternheimer equation
                }
                std::vector<std::complex<double>> rhs = data_.get_dpsi(q_idx, ik, ib);
                if (rhs.empty() || static_cast<int>(dv_sc.size()) != nbands
                    || rhs.size() != dv_sc[ib].size()) {
                    if (dbg) {
                        std::cout << "DBG skip solve ik=" << ik << " ib=" << ib
                                  << " rhs.size=" << rhs.size()
                                  << " dv_sc.size=" << dv_sc.size()
                                  << " dv_sc[ib].size=" << (dv_sc.size() > static_cast<size_t>(ib) ? dv_sc[ib].size() : 999999)
                                  << std::endl;
                    }
                    continue;
                }
                // b = -(dV_ext + dV_sc)|psi_n>
                for (size_t i = 0; i < rhs.size(); ++i) {
                    rhs[i] = -(rhs[i] + dv_sc[ib][i]);
                }
                hamilt_->set_shift(eig_(ik, ib));
                std::vector<std::complex<double>> dpsi_out;
                double res = 0.0;
                stern_.solve(*hamilt_, occ_kq_[ik], rhs, lin_max, lin_thr, dpsi_out, res);
                if (dbg) {
                    double nr = 0.0, nb2 = 0.0;
                    for (size_t i = 0; i < dpsi_out.size(); ++i) {
                        nr += std::norm(dpsi_out[i]);
                        nb2 += std::norm(rhs[i]);
                    }
                    std::cout << "DBG solve ik=" << ik << " ib=" << ib
                              << " eps=" << eig_(ik, ib)
                              << " res=" << res << " |dpsi|=" << std::sqrt(nr)
                              << " |rhs|=" << std::sqrt(nb2)
                              << " finite=" << (std::isfinite(std::sqrt(nr)) ? 1 : 0)
                              << std::endl;
                }
                data_.set_dpsi(q_idx, ik, ib, dpsi_out);
            }
        }

        // ---- 3. first-order density and mixing
        rho_.compute_drho(gs_psi_, wg_, q_idx, data_);
        rho_.mix_drho(q_idx, data_);
        residual = rho_.get_residual(q_idx, data_);
        if (dbg) {
            std::cout << "DBG iter=" << iter << " residual=" << residual
                      << " conv_thr=" << conv_thr_ << std::endl;
        }
        converged = (residual < conv_thr_);
    }
    // stash the converged screened potential and dpsi of this displacement
    // for the two-pass 2n+1 accumulation (term2 cross section needs
    // dV_ext^b + dV_sc^b and dpsi^b of every displacement)
    data_.set_vsc_r(iat, idir, v_sc_r_last);
    {
        std::vector<std::vector<std::vector<std::complex<double>>>> disp(
            nk, std::vector<std::vector<std::complex<double>>>(nbands));
        for (int ik = 0; ik < nk; ++ik) {
            for (int ib = 0; ib < nbands; ++ib) {
                disp[ik][ib] = data_.get_dpsi(q_idx, ik, ib);
            }
        }
        data_.set_dpsi_disp(iat, idir, disp);
    }

    return residual;
}

void DFPT_PW::Impl::solve_pos_resp(int q_idx) {
    // Y^a_{k,v} = P_c x_a|psi_{k,v}> through the Sternheimer equation
    //   (H(k) - eps_v) Y^a_v = P_c [H, x_a]|psi_v>,
    //   [H, x_a]|psi> = -(i/tpiba) dH/dk_a|psi>   (velocity form),
    // exactly the linear solve of QE dvpsi_e (whose rhs negation restores
    // P_c[H,x]psi from commutator_Hx_psi's [x,H] convention). dH/dk_a is the
    // pos_matrix velocity operator: the diagonal kinetic 2 tpiba^2 (k+G)_a
    // plus the separable projector derivative (build_vkb/build_vkb_dk). The
    // solved vector carries the complete conduction-space position response
    // and replaces the empty-eigenvector-truncated r-matrix contraction.
    if (!wired() || hamilt_ == nullptr) {
        return;
    }
    const ModuleBase::Vector3<double> q_cart = data_.get_qvec(q_idx) * ucell_->G;
    const int nk = gs_psi_.get_nk();
    const int nbands = gs_psi_.get_nbands();
    const double tpiba = ucell_->tpiba;
    const double tpiba2 = tpiba * tpiba;
    const int lin_max = data_.get_max_iter();
    const double lin_thr = data_.get_conv_thr();
    const bool dbg = (getenv("DFPT_DEBUG") != nullptr);

    for (int a = 0; a < 3; ++a) {
        std::vector<std::vector<std::vector<std::complex<double>>>> yvec(
            nk, std::vector<std::vector<std::complex<double>>>(nbands));
        for (int ik = 0; ik < nk; ++ik) {
            if (occ_kq_[ik].empty()) {
                continue; // matches the displacement solve guard
            }
            if (last_q_ != q_idx || last_ik_ != ik) {
                hamilt_->set_context(q_cart, ik);
                last_q_ = q_idx;
                last_ik_ = ik;
            }
            const int npwk = pw_wfc_->npwk[ik];
            std::vector<ModuleBase::Vector3<double>> gk(npwk);
            for (int ig = 0; ig < npwk; ++ig) {
                gk[ig] = pw_wfc_->getgpluskcar(ik, ig);
            }
            // dH/dk_a|psi_b> for every band: diagonal kinetic part
            std::vector<std::vector<std::complex<double>>> vel(
                nbands,
                std::vector<std::complex<double>>(npwk, std::complex<double>(0.0, 0.0)));
            for (int ib = 0; ib < nbands; ++ib) {
                for (int ig = 0; ig < npwk; ++ig) {
                    vel[ib][ig] = 2.0 * tpiba2 * gk[ig][a] * gs_psi_(ik, ib, ig);
                }
            }
            // nonlocal derivative part (pos_matrix velocity form; NCPP
            // separable projectors only)
            for (int it = 0; it < ucell_->ntype; ++it) {
                const pseudo& ncpp = ucell_->atoms[it].ncpp;
                const int nh = ncpp.nh;
                if (nh == 0) {
                    continue;
                }
                // projector -> (radial beta index, m channel) table
                std::vector<int> mu_ib(nh, 0);
                std::vector<int> mu_m(nh, 0);
                int mu_idx = 0;
                for (int ib = 0; ib < ncpp.nbeta; ++ib) {
                    const int l = ncpp.lll[ib];
                    for (int m = 0; m < 2 * l + 1; ++m) {
                        if (mu_idx < nh) {
                            mu_ib[mu_idx] = ib;
                            mu_m[mu_idx] = m;
                        }
                        ++mu_idx;
                    }
                }
                for (int ia = 0; ia < ucell_->atoms[it].na; ++ia) {
                    std::vector<std::vector<std::complex<double>>> vkb;
                    pert_.build_vkb(it, ia, gk, vkb);
                    // becp_b[mu] = <vkb_mu|psi_b>
                    std::vector<std::vector<std::complex<double>>> becp(nbands);
                    for (int b = 0; b < nbands; ++b) {
                        becp[b].assign(nh, std::complex<double>(0.0, 0.0));
                        for (int mu = 0; mu < nh; ++mu) {
                            for (int ig = 0; ig < npwk; ++ig) {
                                becp[b][mu] += std::conj(vkb[mu][ig]) * gs_psi_(ik, b, ig);
                            }
                        }
                    }
                    std::vector<std::vector<std::complex<double>>> dvkb;
                    pert_.build_vkb_dk(it, ia, a, gk, vkb, dvkb);
                    // dbecp_b[mu] = <dvkb_mu|psi_b>
                    std::vector<std::vector<std::complex<double>>> dbecp(nbands);
                    for (int b = 0; b < nbands; ++b) {
                        dbecp[b].assign(nh, std::complex<double>(0.0, 0.0));
                        for (int mu = 0; mu < nh; ++mu) {
                            for (int ig = 0; ig < npwk; ++ig) {
                                dbecp[b][mu] += std::conj(dvkb[mu][ig]) * gs_psi_(ik, b, ig);
                            }
                        }
                    }
                    // dV_nl/dk_a|psi_b> = sum_mu |dvkb_mu> (D becp_b)_mu
                    //                              + |vkb_mu> (D dbecp_b)_mu
                    for (int b = 0; b < nbands; ++b) {
                        for (int mu = 0; mu < nh; ++mu) {
                            std::complex<double> out_b(0.0, 0.0);
                            std::complex<double> in_b(0.0, 0.0);
                            for (int nu = 0; nu < nh; ++nu) {
                                if (mu_m[mu] != mu_m[nu]) {
                                    continue;
                                }
                                const double dij = ncpp.dion(mu_ib[mu], mu_ib[nu]);
                                out_b += dij * becp[b][nu];
                                in_b += dij * dbecp[b][nu];
                            }
                            for (int ig = 0; ig < npwk; ++ig) {
                                vel[b][ig] += dvkb[mu][ig] * out_b + vkb[mu][ig] * in_b;
                            }
                        }
                    }
                }
            }
            // solve (H - eps_v) Y = -(i/tpiba) vel for every occupied band
            for (int ib = 0; ib < nbands; ++ib) {
                if (!dfpt_band_occupied(wg_, ik, ib)) {
                    continue;
                }
                std::vector<std::complex<double>> rhs(
                    npwk, std::complex<double>(0.0, 0.0));
                const std::complex<double> fac(0.0, -1.0 / tpiba);
                for (int ig = 0; ig < npwk; ++ig) {
                    rhs[ig] = fac * vel[ib][ig];
                }
                hamilt_->set_shift(eig_(ik, ib));
                double res = 0.0;
                stern_.solve(*hamilt_, occ_kq_[ik], rhs, lin_max, lin_thr,
                             yvec[ik][ib], res);
                if (dbg) {
                    std::cout << "DBG posresp a=" << a << " ik=" << ik
                              << " ib=" << ib << " eps=" << eig_(ik, ib)
                              << " res=" << res << std::endl;
                }
            }
        }
        data_.set_pos_resp(a, yvec);
    }
}

void DFPT_PW::Impl::solve_efield_resp(int q_idx) {
    // E-field SCF response (QE solve_e + dfpt_kernel form): the bare legs
    // Y^a stashed by solve_pos_resp are the field rhs base and the fixed
    // point adds the screened response potential of the mixed drho^E
    // exactly like solve_displacement. The converged dpsi^E,a feeds the
    // SCF dielectric tensor (DFPT_Q0::compute_eps) and the zstar_eu
    // cross-check probe (DFPT_ALEG).
    if (!wired() || hamilt_ == nullptr) {
        return;
    }
    const ModuleBase::Vector3<double> q_cart = data_.get_qvec(q_idx) * ucell_->G;
    const int nrxx = pw_rho_->nrxx;
    const int nk = gs_psi_.get_nk();
    const int nbands = gs_psi_.get_nbands();
    const int lin_max = data_.get_max_iter();
    const double lin_thr = data_.get_conv_thr();

    for (int a = 0; a < 3; ++a) {
        const std::vector<std::vector<std::vector<std::complex<double>>>> yr
            = data_.get_pos_resp(a);
        if (static_cast<int>(yr.size()) != nk) {
            continue; // bare legs not solved: no E response either
        }
        rho_.reset_mixing(q_idx);
        data_.set_drho_g(q_idx, 0,
                         std::vector<std::complex<double>>(pw_rho_->npw,
                                                           std::complex<double>(0.0, 0.0)));
        bool converged = false;
        for (int iter = 0; iter < max_iter_ && !converged; ++iter) {
            // screened response potential of the mixed input density
            // (identical assembly to solve_displacement)
            std::vector<std::complex<double>> v_sc_r(nrxx, std::complex<double>(0.0, 0.0));
            const std::vector<std::complex<double>> drho_in_g = data_.get_drho_g(q_idx, 0);
            if (!drho_in_g.empty() && static_cast<int>(drho_in_g.size()) == pw_rho_->npw) {
                std::vector<std::complex<double>> dv_ha_g;
                rho_.v_hartree_q(q_cart, drho_in_g, dv_ha_g);
                pw_rho_->recip2real(dv_ha_g.data(), v_sc_r.data());
                if (xc_ != nullptr) {
                    std::vector<std::complex<double>> a_r(nrxx);
                    pw_rho_->recip2real(drho_in_g.data(), a_r.data());
                    std::vector<std::complex<double>> b_r;
                    xc_->apply(a_r, b_r);
                    if (static_cast<int>(b_r.size()) == nrxx) {
                        for (int ir = 0; ir < nrxx; ++ir) {
                            v_sc_r[ir] += b_r[ir];
                        }
                    }
                }
            }
            for (int ik = 0; ik < nk; ++ik) {
                if (static_cast<int>(occ_kq_.size()) <= ik || occ_kq_[ik].empty()) {
                    continue;
                }
                std::vector<std::vector<std::complex<double>>> dv_sc;
                pert_.apply_vr(q_idx, ik, v_sc_r, gs_psi_, q_cart, dv_sc);
                if (last_q_ != q_idx || last_ik_ != ik) {
                    hamilt_->set_context(q_cart, ik);
                    last_q_ = q_idx;
                    last_ik_ = ik;
                }
                for (int ib = 0; ib < nbands; ++ib) {
                    if (!dfpt_band_occupied(wg_, ik, ib)) {
                        continue;
                    }
                    if (static_cast<int>(yr[ik][ib].size()) == 0
                        || static_cast<int>(dv_sc.size()) != nbands
                        || yr[ik][ib].size() != dv_sc[ib].size()) {
                        continue;
                    }
                    std::vector<std::complex<double>> rhs(yr[ik][ib].size());
                    for (size_t i = 0; i < rhs.size(); ++i) {
                        rhs[i] = -(yr[ik][ib][i] + dv_sc[ib][i]);
                    }
                    hamilt_->set_shift(eig_(ik, ib));
                    std::vector<std::complex<double>> dpsi_out;
                    double res = 0.0;
                    stern_.solve(*hamilt_, occ_kq_[ik], rhs, lin_max, lin_thr,
                                 dpsi_out, res);
                    data_.set_dpsi(q_idx, ik, ib, dpsi_out);
                }
            }
            rho_.compute_drho(gs_psi_, wg_, q_idx, data_);
            rho_.mix_drho(q_idx, data_);
            const double residual = rho_.get_residual(q_idx, data_);
            converged = (residual < conv_thr_);
            if (converged) {
                std::cout << "DFPT efield dir=" << a
                          << " converged, residual=" << residual
                          << " (iter=" << iter << ")" << std::endl;
            }
        }
        // stash dpsi^E,a before any later solve reuses the slots
        std::vector<std::vector<std::vector<std::complex<double>>>> de(
            nk, std::vector<std::vector<std::complex<double>>>(nbands));
        for (int ik = 0; ik < nk; ++ik) {
            for (int ib = 0; ib < nbands; ++ib) {
                de[ik][ib] = data_.get_dpsi(q_idx, ik, ib);
            }
        }
        data_.set_dpsi_efield(a, de);
    }
}

void DFPT_PW::run() {
    const int nq = pimpl_->qlist_.get_nq();
    for (int q_idx = 0; q_idx < nq; ++q_idx) {
        // Special handling for q=0 (uniform electric field responses):
        // The standard position operator r is ill-defined in periodic systems.
        // Developers should NOT pass a conventional position matrix. Instead,
        // matrix elements should be computed using the well-defined periodic
        // commutator [Ĥ_SCF, r̂]. This is implemented in DFPT_Q0 module.
            if (q_idx == 0 && pimpl_->data_.get_compute_q0()) {
                pimpl_->q0_.compute_q0_response(pimpl_->data_);
            }

        // occupied states at k+q for every k of this q (projector of P_c);
        // also invalidates the shifted-operator context cache
        if (pimpl_->wired()) {
            pimpl_->build_occ_kq(q_idx);
        }

        // position legs of the screened Born charges: the q = 0 Y solves
        // need the projector just built and must land before the two-pass
        // displacement solves below reuse the shifted-operator context
        if (q_idx == 0 && pimpl_->data_.get_compute_q0() && pimpl_->wired()) {
            pimpl_->solve_pos_resp(q_idx);
            // SCF E-field responses of the dielectric tensor: after the
            // bare Y legs they consume, before the displacement solves
            // reuse the slots; the epsilon contraction runs straight after
            // (QE solve_e -> dielec.f90 order)
            pimpl_->solve_efield_resp(q_idx);
            pimpl_->q0_.compute_eps(pimpl_->wg_, pimpl_->data_);
        }

        // Per-irrep self-consistent loop: the little-group irrep
        // decomposition is a placeholder until stage A, so the single
        // available irrep falls back to the full 3N displacement basis.
        // Ledger semantics (B4): one outer pass solves every displacement
        // to its own convergence (solve_displacement restarts each from a
        // zero input density), and the pass residual is the worst final
        // displacement residual; the pass converges when that worst is
        // below conv_thr. An unconverged pass therefore re-runs the full
        // solve, bounded by max_iter_ outer passes, and the residual
        // history keeps an honest record instead of the former
        // unconditional single-pass convergence.
        const int nirr = pimpl_->data_.get_nirr(q_idx);
        for (int irrep = 0; irrep < nirr; ++irrep) {
            pimpl_->data_.set_converged(q_idx, irrep, false);
            pimpl_->data_.set_current_iter(q_idx, irrep, 0);
            while (!pimpl_->data_.get_converged(q_idx, irrep)
                   && pimpl_->data_.get_current_iter(q_idx, irrep) < pimpl_->max_iter_) {
                if (pimpl_->wired()) {
                    const int nat = pimpl_->ucell_->nat;
                    // two passes over the 3N displacement basis: first solve
                    // every displacement to convergence (the 2n+1 accumulation
                    // of displacement b needs the converged dpsi AND screened
                    // potential of every column displacement a), then run the
                    // 2n+1 accumulation for each
                    double worst = 0.0;
                    for (int iat = 0; iat < nat; ++iat) {
                        for (int idir = 0; idir < 3; ++idir) {
                            const double residual = pimpl_->solve_displacement(q_idx, iat, idir);
                            worst = std::max(worst, residual);
                        }
                    }
                    for (int iat = 0; iat < nat; ++iat) {
                        for (int idir = 0; idir < 3; ++idir) {
                            // 2n+1 accumulation of this converged displacement
                            pimpl_->phon_.accumulate_electron(q_idx, iat, idir,
                                                               pimpl_->gs_psi_,
                                                               pimpl_->wg_,
                                                               pimpl_->data_);
                        }
                    }
                    pimpl_->data_.add_residual(q_idx, irrep, worst);
                    pimpl_->data_.set_converged(q_idx, irrep,
                                                worst < pimpl_->data_.get_conv_thr());
                } else {
                    // design-phase skeleton: no bases wired, converge at once
                    pimpl_->data_.add_residual(q_idx, irrep, 0.0);
                    pimpl_->data_.set_converged(q_idx, irrep, true);
                }
                pimpl_->data_.set_current_iter(
                    q_idx, irrep, pimpl_->data_.get_current_iter(q_idx, irrep) + 1);
            }
        }

        // screened Born charges: the Gonze-Lee 2n+1 form consumes the
        // converged (screened) dpsi of every q = 0 displacement stashed by
        // solve_displacement, so it must run after the two-pass solves
        // above and before the LO-TO term below consumes it
        if (q_idx == 0 && pimpl_->data_.get_compute_q0() && pimpl_->wired()) {
            pimpl_->q0_.compute_born(pimpl_->gs_psi_, pimpl_->wg_,
                                     pimpl_->eig_, pimpl_->data_);
        }

        pimpl_->phon_.assemble(q_idx, pimpl_->data_);
        pimpl_->phon_.diagonalize(q_idx, pimpl_->data_);
        if (q_idx == 0 && pimpl_->data_.get_loto()) {
            // non-analytic LO-TO correction along the data-layer direction
            // (default isotropic (1,1,1)/sqrt(3) for cubic crystals;
            // set_loto_dir overrides, e.g. per irrep direction in stage A)
            pimpl_->phon_.add_loto(pimpl_->data_.get_loto_dir(), pimpl_->data_);
            pimpl_->phon_.diagonalize_loto(pimpl_->data_);
        }
    }
}

int DFPT_PW::get_nq() const {
    return pimpl_->qlist_.get_nq();
}

ModuleBase::Vector3<double> DFPT_PW::get_qvec(int q_idx) const {
    return pimpl_->data_.get_qvec(q_idx);
}

std::vector<double> DFPT_PW::get_phonon_freq(int q_idx) const {
    return pimpl_->data_.get_phon_freq(q_idx);
}

std::vector<double> DFPT_PW::get_phon_freq_loto() const {
    return pimpl_->data_.get_phon_freq_loto();
}

ModuleBase::Vector3<double> DFPT_PW::get_loto_dir() const {
    return pimpl_->data_.get_loto_dir();
}

std::string DFPT_PW::format_q_report(int q_idx) const {
    return pimpl_->phon_.format_q_report(q_idx, pimpl_->data_);
}

std::string DFPT_PW::format_loto_report() const {
    return pimpl_->phon_.format_loto_report(pimpl_->data_);
}

ModuleBase::matrix DFPT_PW::get_dielectric_tensor() const {
    return pimpl_->data_.get_dielectric();
}

ModuleBase::matrix DFPT_PW::get_born_charges(int atom_idx) const {
    return pimpl_->data_.get_born(atom_idx);
}

void DFPT_PW::set_qfile(const std::string& filename) {
    pimpl_->qfile_ = filename;
}

void DFPT_PW::set_qmesh(int nqx, int nqy, int nqz) {
    pimpl_->nqx_ = nqx;
    pimpl_->nqy_ = nqy;
    pimpl_->nqz_ = nqz;
}

void DFPT_PW::set_conv_thr(double thr) {
    pimpl_->conv_thr_ = thr;
    pimpl_->data_.set_conv_thr(thr);
}

void DFPT_PW::set_max_iter(int max_iter) {
    pimpl_->max_iter_ = max_iter;
    pimpl_->data_.set_max_iter(max_iter);
}

void DFPT_PW::set_mix_beta(double beta) {
    if (beta > 0.0 && beta <= 1.0) {
        pimpl_->mix_beta_ = beta;
    }
}

void DFPT_PW::set_compute_q0(bool flag) {
    pimpl_->data_.set_compute_q0(flag);
}

void DFPT_PW::set_loto(bool flag) {
    pimpl_->data_.set_loto(flag);
}

void DFPT_PW::set_loto_dir(const ModuleBase::Vector3<double>& dir) {
    pimpl_->data_.set_loto_dir(dir);
}

} // namespace ModuleDFPT
