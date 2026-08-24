// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================

#include "dfpt_pw.h"
#include "dfpt_pw_data.h"
#include "dfpt_irrep_data.h"
#include "dfpt_pert.h"
#include "dfpt_stern.h"
#include "dfpt_rho.h"
#include "dfpt_phon.h"
#include "dfpt_q0.h"
#include "dfpt_metal.h"
#include "dfpt_hamilt_shift.h"
#include "dfpt_kq_basis.h"
#include "source_base/constants.h"
#include <fstream>
#include "source_base/global_function.h"
#include <sstream>
#include "source_cell/qlist.h"
#include "source_pw/module_pwdft/stru_fac.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
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
    ///< BPT debug: empty states at k+q on the k+q G list + their eigenvalues
    std::vector<std::vector<std::vector<std::complex<double>>>> empty_kq_;
    std::vector<std::vector<double>> empty_kq_eig_;
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
    /// with the screening assembly of solve_displacement); design-phase
    /// probe behind DFPT_ALEG for the zstar_eu cross-check
    void solve_efield_resp(int q_idx);

    /// zstar_eu cross-check of the screened Born charges (dpsi^E,scf
    /// contracted with the bare dV^kappa|psi> legs) against the zstar_ue
    /// form compute_born produced (DFPT_ALEG probe)
    void aleg_crosscheck(int q_idx);
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
        // |lambda| ~ 3; the env knob is a design-phase calibration aid
        double mix_beta = pimpl_->mix_beta_;
        if (const char* env_beta = getenv("DFPT_MIX_BETA")) {
            const double parsed = atof(env_beta);
            if (parsed > 0.0 && parsed <= 1.0) {
                mix_beta = parsed;
            }
        }
        pimpl_->rho_.init(nspin, nrxx, pw_rho, pw_wfc, ucell.G, "plain", mix_beta);
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
    // BPT debug companion: empty states at k+q on the same kq ball, for the
    // independent perturbation-theory cross-check of the Sternheimer solve
    // (DFPT_BPT); eig pairs are (eig_(ikq, m), eig_(ik, n))
    empty_kq_.assign(nk, std::vector<std::vector<std::complex<double>>>());
    empty_kq_eig_.assign(nk, std::vector<double>());
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
        const bool want_empty = (getenv("DFPT_BPT") != nullptr);
        int dbg_miss = 0;
        int dbg_tot = 0;
        for (int m = 0; m < nbands; ++m) {
            const bool occ_m = dfpt_band_occupied(wg_, ikq, m);
            if (!occ_m && !want_empty) {
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
                } else {
                    ++dbg_miss;
                }
                ++dbg_tot;
            }
            if (occ_m) {
                occ_kq_[ik].push_back(std::move(state));
            } else {
                empty_kq_[ik].push_back(std::move(state));
                empty_kq_eig_[ik].push_back(eig_(ikq, m));
            }
        }
        if (getenv("DFPT_DEBUG") != nullptr) {
            std::cout << "OCCCHK ik=" << ik << " ikq=" << ikq
                      << " dn=(" << dn_i[0] << "," << dn_i[1] << "," << dn_i[2] << ")"
                      << " npw_kq=" << npw_kq
                      << " npwk_ikq=" << pw_wfc_->npwk[ikq]
                      << " miss=" << dbg_miss << "/" << dbg_tot << std::endl;
            if (dbg_miss > 0 && npw_kq > 0) {
                std::set<std::vector<int>> kq_labels;
                for (int igl = 0; igl < npw_kq; ++igl) {
                    const ModuleBase::Vector3<double> gf = kq.get_gcar(igl) * ginv;
                    kq_labels.insert({static_cast<int>(std::round(gf.x)),
                                      static_cast<int>(std::round(gf.y)),
                                      static_cast<int>(std::round(gf.z))});
                }
                std::set<std::vector<int>> gs_labels;
                std::set<int> gs_igs;
                int ig_max = -1;
                for (int jgl = 0; jgl < pw_wfc_->npwk[ikq]; ++jgl) {
                    const int ig = pw_wfc_->getigl2ig(ikq, jgl);
                    gs_igs.insert(ig);
                    if (ig > ig_max) {
                        ig_max = ig;
                    }
                    const ModuleBase::Vector3<double> gf
                        = pw_wfc_->getgcar(ikq, jgl) * ginv;
                    gs_labels.insert({static_cast<int>(std::round(gf.x)),
                                      static_cast<int>(std::round(gf.y)),
                                      static_cast<int>(std::round(gf.z))});
                }
                std::cout << "OCCCHK   gs_unique_ig=" << gs_igs.size()
                          << "/" << pw_wfc_->npwk[ikq]
                          << " ig_max=" << ig_max
                          << " npw=" << pw_wfc_->npw
                          << " npwk_max=" << pw_wfc_->npwk_max
                          << std::endl;
                int only_kq = 0;
                std::cout << "OCCCHK   labels kq=" << kq_labels.size()
                          << " gs=" << gs_labels.size();
                for (const auto& k : kq_labels) {
                    if (gs_labels.count(k) == 0) {
                        ++only_kq;
                    }
                }
                std::cout << " only_kq=" << only_kq << std::endl;
                int shown = 0;
                for (const auto& k : kq_labels) {
                    if (gs_labels.count(k) == 0) {
                        std::cout << "OCCCHK   kq-only (" << k[0] << "," << k[1] << ","
                                  << k[2] << ")";
                        if (++shown >= 6) {
                            break;
                        }
                    }
                }
                if (shown > 0) {
                    std::cout << std::endl;
                }
                shown = 0;
                for (const auto& k : gs_labels) {
                    if (kq_labels.count(k) == 0) {
                        std::cout << "OCCCHK   gs-only (" << k[0] << "," << k[1] << ","
                                  << k[2] << ")";
                        if (++shown >= 6) {
                            break;
                        }
                    }
                }
                if (shown > 0) {
                    std::cout << std::endl;
                }
                const ModuleBase::Vector3<double> gf0 = kq.get_gcar(0) * ginv;
                std::cout << "OCCCHK   kq key0 gf=(" << gf0.x << "," << gf0.y << ","
                          << gf0.z << ")";
                const ModuleBase::Vector3<double> gj0
                    = pw_wfc_->getgcar(ikq, 0) * ginv;
                std::cout << "  ikq gf0=(" << gj0.x << "," << gj0.y << "," << gj0.z << ")"
                          << std::endl;
            }
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
    // design-phase homogeneous probe: inject a pure A1 trial density on the
    // {200}/{111} shells, drop the external perturbation from the rhs, run
    // one iteration and dump the linear map output M * trial
    const bool jprobe = (getenv("DFPT_JPROBE") != nullptr);
    const bool jprobe_noxc = (getenv("DFPT_JPROBE_NOXC") != nullptr);

    bool converged = false;
    double residual = 0.0;
    const bool dbg = (getenv("DFPT_DEBUG") != nullptr);
    // last screened response potential (hoisted out of the loop: the 2n+1
    // accumulation below needs the converged v_sc of this displacement)
    std::vector<std::complex<double>> v_sc_r_last;
    for (int iter = 0; iter < max_iter_ && !converged; ++iter) {
        data_.set_current_iter(iter);

        if (jprobe && iter == 0) {
            std::vector<std::complex<double>> trial(pw_rho_->npw,
                                                    std::complex<double>(0.0, 0.0));
            double nrm = 0.0;
            for (int ig = 0; ig < pw_rho_->npw; ++ig) {
                const ModuleBase::Vector3<double> gc = pw_rho_->gcar[ig];
                const double g2 = gc * gc;
                const int nax = (std::abs(gc.x) > 1.0e-6)
                                + (std::abs(gc.y) > 1.0e-6)
                                + (std::abs(gc.z) > 1.0e-6);
                if (std::abs(g2 - 4.0) < 1.0e-9 && nax == 1) {
                    trial[ig] = std::complex<double>(1.0, 0.0);
                    nrm += 1.0;
                } else if (std::abs(g2 - 3.0) < 1.0e-9 && nax == 3) {
                    const double sgn = (gc.x * gc.y * gc.z > 0.0) ? 1.0 : -1.0;
                    trial[ig] = std::complex<double>(0.6218, -0.6218 * sgn);
                    nrm += 2.0 * 0.6218 * 0.6218;
                }
            }
            const double inv = 1.0 / std::sqrt(nrm);
            for (size_t i = 0; i < trial.size(); ++i) {
                trial[i] *= inv;
            }
            data_.set_drho_g(q_idx, 0, trial);
        }

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
            const bool noxc = (getenv("DFPT_NOXC") != nullptr);
            double xcs = 1.0;
            if (getenv("DFPT_XCS") != nullptr) {
                xcs = std::atof(getenv("DFPT_XCS"));
            }
            if (xc_ != nullptr && !jprobe_noxc && !noxc) {
                std::vector<std::complex<double>> a_r(nrxx);
                pw_rho_->recip2real(drho_in_g.data(), a_r.data());
                std::vector<std::complex<double>> b_r;
                xc_->apply(a_r, b_r);
                if (static_cast<int>(b_r.size()) == nrxx) {
                    for (int ir = 0; ir < nrxx; ++ir) {
                        v_sc_r[ir] += xcs * b_r[ir];
                    }
                }
                if (getenv("DFPT_MDBG") != nullptr) {
                    static int vdbg_cnt = 0;
                    std::ofstream vf("/tmp/opencode/drho_iters/vsc_"
                                     + std::to_string(vdbg_cnt++) + ".bin",
                                     std::ios::binary);
                    for (int ir = 0; ir < nrxx; ++ir) {
                        vf.write(reinterpret_cast<const char*>(&v_sc_r[ir]),
                                 sizeof(std::complex<double>));
                    }
                    // channel split: rebuild each part for the same input
                    std::vector<std::complex<double>> vh_only(nrxx);
                    pw_rho_->recip2real(dv_ha_g.data(), vh_only.data());
                    std::ofstream hf("/tmp/opencode/drho_iters/vha_"
                                     + std::to_string(vdbg_cnt - 1) + ".bin",
                                     std::ios::binary);
                    for (int ir = 0; ir < nrxx; ++ir) {
                        hf.write(reinterpret_cast<const char*>(&vh_only[ir]),
                                 sizeof(std::complex<double>));
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
                if (getenv("DFPT_MDBG") != nullptr) {
                    static int dump_cnt = 0;
                    std::ofstream df("/tmp/opencode/drho_iters/it"
                                     + std::to_string(dump_cnt++) + ".bin",
                                     std::ios::binary);
                    for (int ig = 0; ig < pw_rho_->npw; ++ig) {
                        df.write(reinterpret_cast<const char*>(&drho_in_g[ig]),
                                 sizeof(std::complex<double>));
                    }
                    static bool dumped_g = false;
                    if (!dumped_g) {
                        dumped_g = true;
                        std::ofstream gf("/tmp/opencode/drho_iters/gcar.bin",
                                         std::ios::binary);
                        for (int ig = 0; ig < pw_rho_->npw; ++ig) {
                            const double gx = pw_rho_->gcar[ig].x;
                            const double gy = pw_rho_->gcar[ig].y;
                            const double gz = pw_rho_->gcar[ig].z;
                            gf.write(reinterpret_cast<const char*>(&gx), sizeof(double));
                            gf.write(reinterpret_cast<const char*>(&gy), sizeof(double));
                            gf.write(reinterpret_cast<const char*>(&gz), sizeof(double));
                        }
                    }
                }
            }
        }
        // DFPT_NOSC: zero the screened part to isolate the bare Sternheimer
        // (design-phase A/B knob for term2 debugging)
        if (getenv("DFPT_NOSC") != nullptr) {
            std::fill(v_sc_r.begin(), v_sc_r.end(), std::complex<double>(0.0, 0.0));
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
                // b = -(dV_ext + dV_sc)|psi_n>; the homogeneous probe keeps
                // only the screening part to isolate the linear map M
                if (jprobe) {
                    for (size_t i = 0; i < rhs.size(); ++i) {
                        rhs[i] = -dv_sc[ib][i];
                    }
                } else {
                    for (size_t i = 0; i < rhs.size(); ++i) {
                        rhs[i] = -(rhs[i] + dv_sc[ib][i]);
                    }
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
                // BPT: independent PT cross-check of the solve for this band.
                // Identity (exact on a complete empty manifold):
                //   <dpsi_n|rhs_n> == sum_m |<psi_m(k+q)|rhs_n>|^2 / (e_m(k+q) - e_n(k))
                // rhs here is the TOTAL (ext+sc) right-hand side with the
                // Sternheimer sign (rhs = -(ext+sc)); both sides use the same
                // object so the sign cancels on the diagonal-identity check.
                if (getenv("DFPT_BPT") != nullptr
                    && static_cast<int>(empty_kq_[ik].size()) > 0
                    && dpsi_out.size() == rhs.size()) {
                    std::complex<double> codedot(0.0, 0.0);
                    double pt = 0.0;
                    double wsum = 0.0;
                    for (size_t i = 0; i < dpsi_out.size(); ++i) {
                        codedot += std::conj(dpsi_out[i]) * rhs[i];
                    }
                    for (size_t im = 0; im < empty_kq_[ik].size(); ++im) {
                        const std::vector<std::complex<double>>& psim
                            = empty_kq_[ik][im];
                        if (psim.size() != rhs.size()) {
                            continue;
                        }
                        std::complex<double> mdot(0.0, 0.0);
                        for (size_t i = 0; i < rhs.size(); ++i) {
                            mdot += std::conj(psim[i]) * rhs[i];
                        }
                        const double denom = empty_kq_eig_[ik][im] - eig_(ik, ib);
                        if (std::abs(denom) > 1.0e-10) {
                            pt += std::norm(mdot) / denom;
                        }
                        wsum += std::norm(mdot);
                    }
                    double nrm2 = 0.0;
                    for (size_t i = 0; i < rhs.size(); ++i) {
                        nrm2 += std::norm(rhs[i]);
                    }
                    std::cout << "BPTCHK q=" << q_idx << " iat=" << iat
                              << " idir=" << idir << " ik=" << ik << " ib=" << ib
                              << " code=(" << codedot.real() << "," << codedot.imag() << ")"
                              << " pt=" << pt
                              << " |<m|rhs>|^2sum=" << wsum
                              << " |rhs|^2=" << nrm2 << std::endl;
                }
            }
        }

        // ---- 3. first-order density and mixing
        rho_.compute_drho(gs_psi_, wg_, q_idx, data_);
        if (jprobe && iter == 0) {
            const std::vector<std::complex<double>> probe_out = data_.get_drho_g(q_idx, 0);
            std::ofstream pf("/tmp/opencode/jprobe_out.bin", std::ios::binary);
            for (int ig = 0; ig < pw_rho_->npw; ++ig) {
                pf.write(reinterpret_cast<const char*>(&probe_out[ig]),
                         sizeof(std::complex<double>));
            }
            std::ofstream vf("/tmp/opencode/jprobe_vsc.bin", std::ios::binary);
            for (int ir = 0; ir < nrxx; ++ir) {
                vf.write(reinterpret_cast<const char*>(&v_sc_r[ir]),
                         sizeof(std::complex<double>));
            }
            pf.close();
            vf.close();
            std::cout << "JPROBE dumped, exiting" << std::endl;
            std::cout.flush();
            std::exit(0);
        }
        rho_.mix_drho(q_idx, data_);
        residual = rho_.get_residual(q_idx, data_);
        data_.add_residual(residual);
        if (dbg) {
            std::cout << "DBG iter=" << iter << " residual=" << residual
                      << " conv_thr=" << conv_thr_ << std::endl;
        }
        converged = (residual < conv_thr_);
    }
    data_.set_converged(converged);
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

    // design-phase validation: dump converged self-consistent drho on the
    // shared real-space grid for direct comparison with finite differences
    if (dbg && q_idx == 0 && iat == 0 && idir == 0) {
        const std::vector<std::complex<double>> dg = data_.get_drho_g(q_idx, 0);
        std::vector<std::complex<double>> dr(pw_rho_->nrxx, std::complex<double>(0.0, 0.0));
        if (static_cast<int>(dg.size()) == pw_rho_->npw) {
            pw_rho_->recip2real(dg.data(), dr.data());
        }
        std::ofstream df("/tmp/opencode/drho_dfpt.dat");
        df << pw_rho_->nx << " " << pw_rho_->ny << " " << pw_rho_->nz << "\n";
        for (int ir = 0; ir < pw_rho_->nrxx; ++ir) {
            df << dr[ir].real() << " " << dr[ir].imag() << "\n";
        }
    }
    // design-phase validation: 8-band perturbation-theory term2 cross-check
    // for the first displacement of the first atom (q = 0, ik = 0 only)
    if (dbg && q_idx == 0 && iat == 0 && idir == 0 && nk > 0 && nbands > 4) {
        // gauge check: <psi_k|dpsi_n> must vanish for occupied k (Sternheimer
        // gauge); a nonzero admixture pollutes term2 via occ-occ dV elements
        for (int ikg = 0; ikg < nk; ++ikg) {
            for (int n = 0; n < 4; ++n) {
                const std::vector<std::complex<double>>& dps = data_.get_dpsi(q_idx, ikg, n);
                const int npwg = gs_psi_.get_nbasis();
                if (static_cast<int>(dps.size()) != npwg) {
                    continue;
                }
                for (int k = 0; k < 4; ++k) {
                    std::complex<double> dot(0.0, 0.0);
                    for (int ig = 0; ig < npwg; ++ig) {
                        dot += std::conj(gs_psi_(ikg, k, ig)) * dps[ig];
                    }
                    std::cout << "PTCHK gauge ik=" << ikg << " n=" << n << " k=" << k
                              << " <psi_k|dpsi_n>=" << dot << std::endl;
                }
            }
        }
        // stash the solved dpsi (apply_dv below reuses the slots)
        std::vector<std::vector<std::complex<double>>> solved(nbands);
        for (int ib = 0; ib < nbands; ++ib) {
            solved[ib] = data_.get_dpsi(q_idx, 0, ib);
        }
        const int npw = gs_psi_.get_nbasis();
        // M[a][m][n] = <psi_m | dV_a | psi_n>
        std::vector<std::vector<std::vector<std::complex<double>>>> mat(
            2, std::vector<std::vector<std::complex<double>>>(
                   nbands, std::vector<std::complex<double>>(nbands, std::complex<double>(0.0, 0.0))));
        for (int a = 0; a < 2; ++a) {
            pert_.build_dv(q_idx, a, idir, data_);
            pert_.apply_dv(q_idx, 0, gs_psi_, data_);
            for (int n = 0; n < nbands; ++n) {
                const std::vector<std::complex<double>> dvpsi = data_.get_dpsi(q_idx, 0, n);
                if (static_cast<int>(dvpsi.size()) != npw) {
                    continue;
                }
                for (int m = 0; m < nbands; ++m) {
                    std::complex<double> dot(0.0, 0.0);
                    for (int ig = 0; ig < npw; ++ig) {
                        dot += std::conj(gs_psi_(0, m, ig)) * dvpsi[ig];
                    }
                    mat[a][m][n] = dot;
                }
            }
        }
        // PT term2 with the 4 empty bands only, b = atom 0 (solved dir)
        for (int a = 0; a < 2; ++a) {
            std::complex<double> pt(0.0, 0.0);
            for (int n = 0; n < 4; ++n) {
                const double w = wg_(0, n);
                if (w < 1.0e-8) {
                    continue;
                }
                for (int m = 4; m < nbands; ++m) {
                    pt += w * std::conj(mat[a][m][n]) * mat[0][m][n]
                          / (eig_(0, n) - eig_(0, m));
                }
            }
            std::cout << "PTCHK term2(a=" << a << ",b=0) PT-4empty=" << 2.0 * pt << std::endl;
        }
        // element dump: M[a][m][n] for m empty, n occupied
        for (int a = 0; a < 2; ++a) {
            for (int m = 4; m < nbands; ++m) {
                for (int n = 0; n < 4; ++n) {
                    std::cout << "PTCHK M a=" << a << " m=" << m << " n=" << n
                              << " (" << mat[a][m][n].real() << "," << mat[a][m][n].imag() << ")"
                              << std::endl;
                }
            }
        }
        // solved-dpsi band projections <psi_m|dpsi_n>
        for (int n = 0; n < 4; ++n) {
            for (int m = 0; m < nbands; ++m) {
                std::complex<double> dot(0.0, 0.0);
                for (int ig = 0; ig < npw; ++ig) {
                    dot += std::conj(gs_psi_(0, m, ig)) * solved[n][ig];
                }
                std::cout << "PTCHK proj n=" << n << " m=" << m
                          << " <psi_m|dpsi_n>=(" << dot.real() << "," << dot.imag() << ")"
                          << std::endl;
            }
        }
        // Hellmann-Feynman check targets: <psi_n|dV_a|psi_n> vs FD of eps_n
        for (int a = 0; a < 2; ++a) {
            for (int n = 0; n < 4; ++n) {
                std::cout << "PTCHK HF a=" << a << " n=" << n
                          << " <psi|dV|psi>=" << mat[a][n][n] << std::endl;
            }
        }
        // design-phase validation: rebuild the converged screened potential
        // and compare de_code(n) = <dV_ext> + <dv_sc> against FD eigenvalue
        // derivatives (localizes response errors in the screening channel)
        {
            const std::vector<std::complex<double>> dg_in = data_.get_drho_g(q_idx, 0);
            if (static_cast<int>(dg_in.size()) == pw_rho_->npw) {
                std::vector<std::complex<double>> vha_g;
                rho_.v_hartree_q(q_cart, dg_in, vha_g);
                std::vector<std::complex<double>> v_sc_r2(pw_rho_->nrxx, std::complex<double>(0.0, 0.0));
                std::vector<std::complex<double>> vh_r(pw_rho_->nrxx);
                pw_rho_->recip2real(vha_g.data(), vh_r.data());
                std::vector<std::complex<double>> vx_r;
            if (xc_ != nullptr) {
                std::vector<std::complex<double>> ar(pw_rho_->nrxx);
                    pw_rho_->recip2real(dg_in.data(), ar.data());
                    xc_->apply(ar, vx_r);
                }
                for (int ir = 0; ir < pw_rho_->nrxx; ++ir) {
                    v_sc_r2[ir] = vh_r[ir];
                    if (static_cast<int>(vx_r.size()) == pw_rho_->nrxx) {
                        v_sc_r2[ir] += vx_r[ir];
                    }
                }
                std::vector<std::vector<std::complex<double>>> dvsc;
                pert_.apply_vr(q_idx, 0, v_sc_r2, gs_psi_, q_cart, dvsc);
                for (int n = 0; n < nbands; ++n) {
                    const std::vector<std::complex<double>>& v = dvsc[n];
                    if (static_cast<int>(v.size()) != npw) {
                        continue;
                    }
                    std::complex<double> dot(0.0, 0.0);
                    for (int ig = 0; ig < npw; ++ig) {
                        dot += std::conj(gs_psi_(0, n, ig)) * v[ig];
                    }
                    std::cout << "PTCHK de n=" << n
                              << " <psi|dv_sc|psi>=(" << dot.real() << "," << dot.imag() << ")"
                              << " de_code=" << (mat[0][n][n] + dot).real() << std::endl;
                }
            }
        }
        // restore the solved dpsi
        for (int ib = 0; ib < nbands; ++ib) {
            if (!solved[ib].empty()) {
                data_.set_dpsi(q_idx, 0, ib, solved[ib]);
            }
        }
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
                    // DKCHK: analytic dk-derivative vs central difference of
                    // the GS-validated build_vkb (first atom with nh > 0 only)
                    if (getenv("DFPT_DKCHK") != nullptr && ia == 0 && ik == 0) {
                        const double dd = 1.0e-5;
                        std::vector<ModuleBase::Vector3<double>> gk_p(npwk);
                        std::vector<ModuleBase::Vector3<double>> gk_m(npwk);
                        for (int ig = 0; ig < npwk; ++ig) {
                            gk_p[ig] = gk[ig];
                            gk_m[ig] = gk[ig];
                            gk_p[ig][a] += dd;
                            gk_m[ig][a] -= dd;
                        }
                        std::vector<std::vector<std::complex<double>>> vkb_p;
                        std::vector<std::vector<std::complex<double>>> vkb_m;
                        pert_.build_vkb(it, ia, gk_p, vkb_p);
                        pert_.build_vkb(it, ia, gk_m, vkb_m);
                        for (int mu = 0; mu < nh; ++mu) {
                            std::complex<double> ddot(0.0, 0.0);
                            double ndk = 0.0;
                            double nnum = 0.0;
                            for (int ig = 0; ig < npwk; ++ig) {
                                const std::complex<double> num
                                    = (vkb_p[mu][ig] - vkb_m[mu][ig]) / (2.0 * dd);
                                ddot += std::conj(dvkb[mu][ig]) * num;
                                ndk += std::norm(dvkb[mu][ig]);
                                nnum += std::norm(num);
                            }
                            std::cout << "DKCHK it=" << it << " a=" << a
                                      << " mu=" << mu
                                      << " <dk|num>=" << ddot
                                      << " |dk|^2=" << ndk
                                      << " |num|^2=" << nnum << std::endl;
                        }
                    }
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
            const bool ychk = (getenv("DFPT_YCHK") != nullptr);
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
                if (ychk && npwk > 0) {
                    // eigen-projection identity: <u_m|Y_v> must equal the
                    // velocity-form matrix element -i<u_m|dH/dk_a|v> /
                    // (tpiba (eps_m - eps_v)) for every nondegenerate
                    // conduction m (pos_matrix cross-check)
                    for (int m = 0; m < nbands; ++m) {
                        if (dfpt_band_occupied(wg_, ik, m)) {
                            continue;
                        }
                        const double de = eig_(ik, m) - eig_(ik, ib);
                        if (std::abs(de) < 1.0e-8) {
                            continue;
                        }
                        std::complex<double> pdot(0.0, 0.0);
                        std::complex<double> ydot(0.0, 0.0);
                        for (int ig = 0; ig < npwk; ++ig) {
                            pdot += std::conj(gs_psi_(ik, m, ig)) * vel[ib][ig];
                            if (ig < static_cast<int>(yvec[ik][ib].size())) {
                                ydot += std::conj(gs_psi_(ik, m, ig))
                                        * yvec[ik][ib][ig];
                            }
                        }
                        const std::complex<double> expect
                            = std::complex<double>(0.0, -1.0) * pdot
                              / (tpiba * de);
                        std::cout << "YCHK a=" << a << " ik=" << ik
                                  << " v=" << ib << " m=" << m
                                  << " <m|Y>=(" << ydot.real() << "," << ydot.imag() << ")"
                                  << " X_mv=(" << expect.real() << "," << expect.imag() << ")"
                                  << " res=" << res << std::endl;
                    }
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
            data_.set_current_iter(iter);
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

void DFPT_PW::Impl::aleg_crosscheck(int q_idx) {
    // zstar_eu cross-check (QE zstar_eu.f90):
    //   Z*(E,Us)_iat(a, d) = zion delta_ad
    //     - 2 sum_k w_k Re sum_v <dpsi^E,a(scf) | dV^kappa(iat,d)_bare | psi_v>
    // with the bare kappa leg from the same build_dv/apply_dv pair the
    // displacement solves consumed. Printed against the plain (no star
    // rotation) zstar_ue form from the dpsi_disp stash and the stored
    // star-rotated Born charges, plus the SCF dielectric tensor
    //   eps = 1 - (16 pi / omega) sum_k wg Re <Y^a | dpsi^E,b>
    // which is the same contraction compute_eps produces in production.
    // A == B incriminates
    // a shared operand; A clean incriminates the compute_born contraction.
    const int nat = ucell_->nat;
    const int nk = gs_psi_.get_nk();
    const int nbands = gs_psi_.get_nbands();
    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> de(3);
    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> yr(3);
    for (int a = 0; a < 3; ++a) {
        de[a] = data_.get_dpsi_efield(a);
        yr[a] = data_.get_pos_resp(a);
        if (static_cast<int>(de[a].size()) != nk
            || static_cast<int>(yr[a].size()) != nk) {
            std::cout << "ALEG: missing E/position response stashes" << std::endl;
            return;
        }
    }
    // A-leg: plain k sum with the full wg weight (QE form). Complex
    // accumulation: a spurious relative phase between the E and kappa legs
    // cancels in the EE and kappa-kappa channels but rotates this cross
    // channel, so Im(z) vs Re(z) is the phase detector (ASR: |z| = zion/2)
    std::vector<ModuleBase::matrix> za(nat, ModuleBase::matrix(3, 3, true));
    std::vector<std::vector<std::complex<double>>> za_c(
        nat, std::vector<std::complex<double>>(9, std::complex<double>(0.0, 0.0)));
    for (int iat = 0; iat < nat; ++iat) {
        for (int idir = 0; idir < 3; ++idir) {
            pert_.build_dv(q_idx, iat, idir, data_);
            for (int ik = 0; ik < nk; ++ik) {
                if (static_cast<int>(occ_kq_.size()) <= ik || occ_kq_[ik].empty()) {
                    continue;
                }
                pert_.apply_dv(q_idx, ik, gs_psi_, data_);
                for (int v = 0; v < nbands; ++v) {
                    if (!dfpt_band_occupied(wg_, ik, v)) {
                        continue;
                    }
                    const std::vector<std::complex<double>> b
                        = data_.get_dpsi(q_idx, ik, v);
                    const int npw = static_cast<int>(b.size());
                    if (npw <= 0) {
                        continue;
                    }
                    for (int a = 0; a < 3; ++a) {
                        if (static_cast<int>(de[a][ik][v].size()) != npw) {
                            continue;
                        }
                        std::complex<double> dot(0.0, 0.0);
                        for (int ig = 0; ig < npw; ++ig) {
                            dot += std::conj(de[a][ik][v][ig]) * b[ig];
                        }
                        za[iat](a, idir) += wg_(ik, v) * dot.real();
                        za_c[iat][3 * a + idir] += wg_(ik, v) * dot;
                        if (a == 0 && idir == 0 && ik < 3) {
                            std::cout << "ALEG-K A ik=" << ik << " v=" << v
                                      << " z=(" << dot.real() << "," << dot.imag() << ")"
                                      << " wg=" << wg_(ik, v) << std::endl;
                        }
                    }
                }
            }
        }
    }
    // plain B-leg from the displacement stashes (same sum, no stars)
    std::vector<ModuleBase::matrix> zb(nat, ModuleBase::matrix(3, 3, true));
    for (int iat = 0; iat < nat; ++iat) {
        for (int idir = 0; idir < 3; ++idir) {
            const std::vector<std::vector<std::vector<std::complex<double>>>> disp
                = data_.get_dpsi_disp(iat, idir);
            if (static_cast<int>(disp.size()) != nk) {
                continue;
            }
            for (int ik = 0; ik < nk; ++ik) {
                for (int v = 0; v < nbands; ++v) {
                    if (!dfpt_band_occupied(wg_, ik, v)) {
                        continue;
                    }
                    const int npw = static_cast<int>(disp[ik][v].size());
                    if (npw <= 0) {
                        continue;
                    }
                    for (int a = 0; a < 3; ++a) {
                        if (static_cast<int>(yr[a][ik][v].size()) != npw) {
                            continue;
                        }
                        std::complex<double> dot(0.0, 0.0);
                        for (int ig = 0; ig < npw; ++ig) {
                            dot += std::conj(disp[ik][v][ig]) * yr[a][ik][v][ig];
                        }
                        zb[iat](a, idir) += wg_(ik, v) * dot.real();
                        if (a == 0 && idir == 0 && ik < 3) {
                            std::cout << "ALEG-K B ik=" << ik << " v=" << v
                                      << " z=(" << dot.real() << "," << dot.imag() << ")"
                                      << " wg=" << wg_(ik, v) << std::endl;
                        }
                    }
                }
            }
        }
    }
    // SCF dielectric tensor (complex accumulation: Im(chi) = phase detector)
    ModuleBase::matrix eps(3, 3, true);
    std::vector<std::complex<double>> eps_c(9, std::complex<double>(0.0, 0.0));
    for (int ik = 0; ik < nk; ++ik) {
        for (int v = 0; v < nbands; ++v) {
            if (!dfpt_band_occupied(wg_, ik, v)) {
                continue;
            }
            for (int a = 0; a < 3; ++a) {
                const int npw = static_cast<int>(yr[a][ik][v].size());
                if (npw <= 0) {
                    continue;
                }
                for (int b = 0; b < 3; ++b) {
                    if (static_cast<int>(de[b][ik][v].size()) != npw) {
                        continue;
                    }
                    std::complex<double> dot(0.0, 0.0);
                    for (int ig = 0; ig < npw; ++ig) {
                        dot += std::conj(yr[a][ik][v][ig]) * de[b][ik][v][ig];
                    }
                    eps(a, b) += wg_(ik, v) * dot.real();
                    eps_c[3 * a + b] += wg_(ik, v) * dot;
                }
            }
        }
    }
    // PTCROSS: bare cross-form operator diagnostic. A fresh bare E-leg
    // solve x = M^-1(-Y^0) (rhs fully known, no screening) is contracted
    // with the bare kappa leg b^kappa: <x|b> must equal the spectral sum
    // over the available empty bands -<Y|m><m|b>/(eps_m-eps_v); the
    // mismatch beyond that band space measures what the truncated
    // manifold misses in the cross channel that diagonal norm checks
    // cannot see. The kappa-side solve x^kappa = M^-1(-b) gives the
    // hermiticity mirror <x^kappa|Y> = conj(<x|Y^0-weighted b>).
    if (getenv("DFPT_PTCROSS") != nullptr) {
        const ModuleBase::Vector3<double> q_cart_p
            = data_.get_qvec(q_idx) * ucell_->G;
        const int ikmax = std::min(nk, 3);
        for (int idir = 0; idir < 3; ++idir) {
            pert_.build_dv(q_idx, 0, idir, data_);
            for (int ik = 0; ik < ikmax; ++ik) {
                if (static_cast<int>(occ_kq_.size()) <= ik || occ_kq_[ik].empty()) {
                    continue;
                }
                pert_.apply_dv(q_idx, ik, gs_psi_, data_);
                if (last_q_ != q_idx || last_ik_ != ik) {
                    hamilt_->set_context(q_cart_p, ik);
                    last_q_ = q_idx;
                    last_ik_ = ik;
                }
                for (int v = 0; v < nbands; ++v) {
                    if (!dfpt_band_occupied(wg_, ik, v)) {
                        continue;
                    }
                    const int npw = static_cast<int>(yr[0][ik][v].size());
                    if (npw <= 0) {
                        continue;
                    }
                    const std::vector<std::complex<double>> b
                        = data_.get_dpsi(q_idx, ik, v);
                    if (static_cast<int>(b.size()) != npw) {
                        continue;
                    }
                    // empty-band overlaps at this k
                    std::vector<std::complex<double>> myv;
                    std::vector<std::complex<double>> mbv;
                    std::vector<double> dev;
                    double wsum = 0.0;
                    for (int m = 0; m < nbands; ++m) {
                        if (dfpt_band_occupied(wg_, ik, m)) {
                            continue;
                        }
                        const double de = eig_(ik, m) - eig_(ik, v);
                        if (std::abs(de) < 1.0e-8) {
                            continue;
                        }
                        std::complex<double> my(0.0, 0.0);
                        std::complex<double> mb(0.0, 0.0);
                        for (int i = 0; i < npw; ++i) {
                            my += std::conj(gs_psi_(ik, m, i)) * yr[0][ik][v][i];
                            mb += std::conj(gs_psi_(ik, m, i)) * b[i];
                        }
                        myv.push_back(my);
                        mbv.push_back(mb);
                        dev.push_back(de);
                        wsum += std::norm(my);
                    }
                    // E-side solve: (H - eps_v) x = -Y^0_v
                    std::vector<std::complex<double>> rhsE(npw);
                    for (int i = 0; i < npw; ++i) {
                        rhsE[i] = -yr[0][ik][v][i];
                    }
                    hamilt_->set_shift(eig_(ik, v));
                    std::vector<std::complex<double>> xE;
                    double resE = 0.0;
                    stern_.solve(*hamilt_, occ_kq_[ik], rhsE,
                                 data_.get_max_iter(), data_.get_conv_thr(),
                                 xE, resE);
                    // kappa-side solve: (H - eps_v) xk = -b
                    std::vector<std::complex<double>> rhsK(npw);
                    for (int i = 0; i < npw; ++i) {
                        rhsK[i] = -b[i];
                    }
                    hamilt_->set_shift(eig_(ik, v));
                    std::vector<std::complex<double>> xK;
                    double resK = 0.0;
                    stern_.solve(*hamilt_, occ_kq_[ik], rhsK,
                                 data_.get_max_iter(), data_.get_conv_thr(),
                                 xK, resK);
                    std::complex<double> crossE(0.0, 0.0);
                    std::complex<double> crossK(0.0, 0.0);
                    std::complex<double> diagE(0.0, 0.0);
                    double pt_cross = 0.0;
                    double pt_diagE = 0.0;
                    double pt_diagK = 0.0;
                    for (int i = 0; i < npw; ++i) {
                        if (static_cast<int>(xE.size()) == npw) {
                            crossE += std::conj(xE[i]) * b[i];
                            diagE += std::conj(xE[i]) * yr[0][ik][v][i];
                        }
                        if (static_cast<int>(xK.size()) == npw) {
                            crossK += std::conj(xK[i]) * yr[0][ik][v][i];
                        }
                    }
                    for (size_t im = 0; im < myv.size(); ++im) {
                        pt_cross += (-std::conj(myv[im]) * mbv[im]
                                     / dev[im]).real();
                        pt_diagE += -std::norm(myv[im]) / dev[im];
                        pt_diagK += -std::norm(mbv[im]) / dev[im];
                    }
                    std::cout << "PTCROSS d=" << idir << " ik=" << ik
                              << " v=" << v
                              << " crossE=(" << crossE.real() << ","
                              << crossE.imag() << ")"
                              << " crossK=(" << crossK.real() << ","
                              << crossK.imag() << ")"
                              << " pt=" << pt_cross
                              << " rE=" << (pt_cross != 0.0
                                                    ? crossE.real() / pt_cross
                                                    : 0.0)
                              << " diagE=(" << diagE.real() << ","
                              << diagE.imag() << ") ptDiagE=" << pt_diagE
                              << " ptDiagK=" << pt_diagK
                              << " wsumY=" << wsum
                              << " resE=" << resE << " resK=" << resK
                              << std::endl;
                }
            }
        }
    }
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            eps(a, b) *= -16.0 * ModuleBase::PI / ucell_->omega;
            if (a == b) {
                eps(a, b) += 1.0;
            }
        }
    }
    std::cout << "ALEG eps_scf:" << std::endl;
    for (int a = 0; a < 3; ++a) {
        std::cout << "  " << eps(a, 0) << " " << eps(a, 1) << " " << eps(a, 2)
                  << std::endl;
    }
    std::cout << "ALEG chi_scf complex diag: (16pi/omega)*z per component"
              << std::endl;
    for (int a = 0; a < 3; ++a) {
        const std::complex<double> z = eps_c[4 * a] * (-16.0 * ModuleBase::PI
                                                       / ucell_->omega);
        std::cout << "  a=" << a << " z=(" << z.real() << "," << z.imag()
                  << ") |z|=" << std::abs(z) << std::endl;
    }
    const ModuleBase::matrix eps_ipa = data_.get_dielectric();
    std::cout << "ALEG eps_ipa(stored):" << std::endl;
    for (int a = 0; a < 3; ++a) {
        std::cout << "  " << eps_ipa(a, 0) << " " << eps_ipa(a, 1) << " "
                  << eps_ipa(a, 2) << std::endl;
    }
    for (int iat = 0; iat < nat; ++iat) {
        const int it = ucell_->iat2it[iat];
        const double zion = ucell_->atoms[it].ncpp.zv;
        std::cout << "ALEG atom " << iat << " zion=" << zion << std::endl;
        std::cout << "  A(eu plain):" << std::endl;
        for (int a = 0; a < 3; ++a) {
            std::cout << "    ";
            for (int d = 0; d < 3; ++d) {
                std::cout << ((a == d) ? zion : 0.0) - 2.0 * za[iat](a, d) << " ";
            }
            std::cout << std::endl;
        }
        std::cout << "  A complex diag: z00..z22 (sum wg<de|b>)" << std::endl;
        for (int a = 0; a < 3; ++a) {
            const std::complex<double> z = za_c[iat][4 * a];
            std::cout << "    z" << a << a << " = (" << z.real() << ","
                      << z.imag() << ") |z|=" << std::abs(z)
                      << " |z|/(zion/2)=" << std::abs(z) / (0.5 * zion)
                      << std::endl;
        }
        std::cout << "  B(ue plain):" << std::endl;
        for (int a = 0; a < 3; ++a) {
            std::cout << "    ";
            for (int d = 0; d < 3; ++d) {
                std::cout << ((a == d) ? zion : 0.0) - 2.0 * zb[iat](a, d) << " ";
            }
            std::cout << std::endl;
        }
        const ModuleBase::matrix zstar = data_.get_born(iat);
        std::cout << "  B(ue star-rot stored):" << std::endl;
        for (int a = 0; a < 3; ++a) {
            std::cout << "    " << zstar(a, 0) << " " << zstar(a, 1) << " "
                      << zstar(a, 2) << std::endl;
        }
    }
}

void DFPT_PW::run() {
    const int nq = pimpl_->qlist_.get_nq();
    DFPT_IrrepData irrep_data(pimpl_->data_);
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
        const int nirr = irrep_data.get_nirr(q_idx);
        for (int irrep = 0; irrep < nirr; ++irrep) {
            irrep_data.set_converged(q_idx, irrep, false);
            irrep_data.set_current_iter(q_idx, irrep, 0);
            while (!irrep_data.get_converged(q_idx, irrep)
                   && irrep_data.get_current_iter(q_idx, irrep) < pimpl_->max_iter_) {
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
                    irrep_data.add_residual(q_idx, irrep, worst);
                } else {
                    // design-phase skeleton: no bases wired, converge at once
                    irrep_data.add_residual(q_idx, irrep, 0.0);
                }
                irrep_data.set_converged(q_idx, irrep, true);
            }
        }

        // screened Born charges: the Gonze-Lee 2n+1 form consumes the
        // converged (screened) dpsi of every q = 0 displacement stashed by
        // solve_displacement, so it must run after the two-pass solves
        // above and before the LO-TO term below consumes it
        if (q_idx == 0 && pimpl_->data_.get_compute_q0() && pimpl_->wired()) {
            pimpl_->q0_.compute_born(pimpl_->gs_psi_, pimpl_->wg_,
                                     pimpl_->eig_, pimpl_->data_);
            if (getenv("DFPT_ALEG") != nullptr) {
                pimpl_->aleg_crosscheck(q_idx);
            }
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
