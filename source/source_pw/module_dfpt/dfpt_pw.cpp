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
#include "source_base/global_function.h"
#include "source_cell/qlist.h"
#include "source_pw/module_pwdft/stru_fac.h"

#include <algorithm>
#include <cmath>
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

    int nqx_ = 1, nqy_ = 1, nqz_ = 1;
    double conv_thr_ = 1e-8;
    int max_iter_ = 100;

    bool wired() const { return pw_rho_ != nullptr && pw_wfc_ != nullptr; }

    /// occupied-state projector set at k+q for every k of this q (commensurate
    /// q: kvec_d[ik] + q must be a k point of the ground-state list mod lattice)
    void build_occ_kq(int q_idx);

    /// one self-consistent Sternheimer cycle for the displacement (iat, idir)
    /// at q; returns the achieved density residual (zero when unwired)
    double solve_displacement(int q_idx, int iat, int idir);
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
    pimpl_->xc_ = xc;
    pimpl_->nelec_ = nelec;
    pimpl_->ecutwfc_ = ecutwfc;
    pimpl_->dftu_ = dftu;

    std::vector<int> mp_grid = {pimpl_->nqx_, pimpl_->nqy_, pimpl_->nqz_};
    pimpl_->qlist_.generate_mesh(ucell, ucell.symm, mp_grid, true);

    int nq = pimpl_->qlist_.get_nq();
    int nk = psi.get_nk();
    int nbands = psi.get_nbands();
    int npw_max = psi.get_current_ngk();
    int nrxx = (pw_rho != nullptr) ? pw_rho->nrxx : 0;
    int nspin = 1;
    int nat = ucell.nat;

    if (pw_rho != nullptr && pw_wfc != nullptr && sf != nullptr) {
        pimpl_->pert_.init(ucell, pw_rho, pw_wfc, *sf);
        pimpl_->rho_.init(nspin, nrxx, pw_rho, pw_wfc, ucell.G, "plain", 0.7);
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
            ModuleBase::WARNING_QUIT("DFPT_PW::build_occ_kq",
                                     "k+q is not a point of the ground-state k list: "
                                     "the DFPT q mesh must be commensurate with the "
                                     "k mesh (and inside the first Brillouin zone).");
        }

        DFPT_KQ_Basis kq;
        kq.init(pw_wfc_, q_cart, ik);
        const int npw_kq = kq.get_npwk();

        // reverse map FFT cell -> per-k G index at ikq (the two k balls are
        // enumerated per k, only the cell position identifies the G)
        std::vector<int> jgl_of_cell(pw_wfc_->nxyz, -1);
        for (int jgl = 0; jgl < pw_wfc_->npwk[ikq]; ++jgl) {
            const int isz = pw_wfc_->getigl2isz(ikq, jgl);
            const int iz = isz % pw_wfc_->nz;
            const int is = isz / pw_wfc_->nz;
            const int ixy = pw_wfc_->is2fftixy[is];
            const int ix = ixy / pw_wfc_->fftny;
            const int iy = ixy % pw_wfc_->fftny;
            jgl_of_cell[(ix * pw_wfc_->ny + iy) * pw_wfc_->nz + iz] = jgl;
        }

        const int nbands = gs_psi_.get_nbands();
        for (int m = 0; m < nbands; ++m) {
            if (wg_(ikq, m) < 1.0e-8) {
                continue; // empty at k+q: outside the P_c projector
            }
            std::vector<std::complex<double>> state(npw_kq, std::complex<double>(0.0, 0.0));
            for (int igl = 0; igl < npw_kq; ++igl) {
                const int isz = kq.get_ig2isz(igl);
                const int iz = isz % pw_wfc_->nz;
                const int is = isz / pw_wfc_->nz;
                const int ixy = pw_wfc_->is2fftixy[is];
                const int ix = ixy / pw_wfc_->fftny;
                const int iy = ixy % pw_wfc_->fftny;
                const int jgl = jgl_of_cell[(ix * pw_wfc_->ny + iy) * pw_wfc_->nz + iz];
                if (jgl >= 0) {
                    state[igl] = gs_psi_(ikq, m, jgl);
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

    const int lin_max = data_.get_max_iter();
    const double lin_thr = data_.get_conv_thr();

    bool converged = false;
    double residual = 0.0;
    for (int iter = 0; iter < max_iter_ && !converged; ++iter) {
        data_.set_current_iter(iter);

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
        }

        // ---- 2. Sternheimer solve of every occupied (k, band)
        for (int ik = 0; ik < nk; ++ik) {
            if (static_cast<int>(occ_kq_.size()) <= ik || occ_kq_[ik].empty()) {
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
            }
            for (int ib = 0; ib < nbands; ++ib) {
                if (wg_(ik, ib) < 1.0e-8) {
                    continue; // unoccupied: no Sternheimer equation
                }
                std::vector<std::complex<double>> rhs = data_.get_dpsi(q_idx, ik, ib);
                if (rhs.empty() || static_cast<int>(dv_sc.size()) != nbands
                    || rhs.size() != dv_sc[ib].size()) {
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
                data_.set_dpsi(q_idx, ik, ib, dpsi_out);
            }
        }

        // ---- 3. first-order density and mixing
        rho_.compute_drho(gs_psi_, wg_, q_idx, data_);
        rho_.mix_drho(q_idx, data_);
        residual = rho_.get_residual(q_idx, data_);
        data_.add_residual(residual);
        converged = (residual < conv_thr_);
    }
    data_.set_converged(converged);
    return residual;
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
            if (pimpl_->wired()) {
                // the dielectric tensor and Born charges of the C6 velocity
                // form; the LO-TO term below consumes them
                pimpl_->q0_.compute_eps(pimpl_->gs_psi_, pimpl_->wg_, pimpl_->eig_, pimpl_->data_);
                pimpl_->q0_.compute_born(pimpl_->gs_psi_, pimpl_->wg_, pimpl_->eig_, pimpl_->data_);
            }
        }

        // occupied states at k+q for every k of this q (projector of P_c);
        // also invalidates the shifted-operator context cache
        if (pimpl_->wired()) {
            pimpl_->build_occ_kq(q_idx);
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
                    double worst = 0.0;
                    for (int iat = 0; iat < nat; ++iat) {
                        for (int idir = 0; idir < 3; ++idir) {
                            const double residual = pimpl_->solve_displacement(q_idx, iat, idir);
                            worst = std::max(worst, residual);
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

        pimpl_->phon_.assemble(q_idx, pimpl_->data_);
        pimpl_->phon_.diagonalize(q_idx, pimpl_->data_);
        if (q_idx == 0 && pimpl_->data_.get_loto()) {
            // non-analytic LO-TO correction along a documented default
            // direction (isotropic for cubic crystals; a general q->0
            // direction control arrives with the irrep machinery of stage A)
            const double inv = 1.0 / std::sqrt(3.0);
            pimpl_->phon_.add_loto(ModuleBase::Vector3<double>(inv, inv, inv), pimpl_->data_);
        }
    }
}

std::vector<double> DFPT_PW::get_phonon_freq(int q_idx) const {
    return pimpl_->data_.get_phon_freq(q_idx);
}

ModuleBase::matrix DFPT_PW::get_dielectric_tensor() const {
    return pimpl_->data_.get_dielectric();
}

ModuleBase::matrix DFPT_PW::get_born_charges(int atom_idx) const {
    return pimpl_->data_.get_born(atom_idx);
}

void DFPT_PW::set_parameters(const std::string& param_file) {
    (void)param_file;
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

} // namespace ModuleDFPT
