// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in the design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================

#include "dfpt_rho.h"

#include "dfpt_kq_basis.h"
#include "source_base/constants.h"
#include "source_base/global_function.h"
#include "source_base/module_mixing/plain_mixing.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <functional>
#include <vector>

namespace ModuleDFPT {

DFPT_Rho::DFPT_Rho() {}

DFPT_Rho::~DFPT_Rho() {
    if (mixer_ != nullptr) {
        delete mixer_;
        mixer_ = nullptr;
    }
}

void DFPT_Rho::init(int nspin, int nrxx, ModulePW::PW_Basis* pw_rho,
                    ModulePW::PW_Basis_K* pw_wfc,
                    const ModuleBase::Matrix3& recip_matrix,
                    const std::string& mix_type, double mix_beta) {
    nspin_ = nspin;
    nrxx_ = nrxx;
    pw_rho_ = pw_rho;
    pw_wfc_ = pw_wfc;
    recip_matrix_ = recip_matrix;
    mix_beta_ = mix_beta;
    if (mix_type != "plain")
    {
        ModuleBase::WARNING_QUIT("DFPT_Rho",
                                 "only plain mixing is supported in the design phase");
    }
    delete mixer_;
    mixer_ = new Base_Mixing::Plain_Mixing(mix_beta_);
}

void DFPT_Rho::compute_drho(const psi::Psi<std::complex<double>>& psi,
                            const ModuleBase::matrix& wg, int q_idx,
                            DFPT_PW_Data& data) {
    if (pw_rho_ == nullptr || pw_wfc_ == nullptr) {
        return;
    }
    if (nspin_ != 1)
    {
        ModuleBase::WARNING_QUIT("DFPT_Rho",
                                 "only nspin = 1 is supported in the design phase");
    }
    const int nk = psi.get_nk();
    const int nbands = psi.get_nbands();
    const ModuleBase::Vector3<double> q_frac = data.get_qvec(q_idx);
    const ModuleBase::Vector3<double> q_cart = q_frac * recip_matrix_;

    // rho-ig -> FFT-cell reverse map through the shared (ix,iy,iz) triple
    // (the rho/wfc stick encodings are not interchangeable, C1 finding)
    std::vector<int> ig_of_cell(pw_rho_->nxyz, -1);
    for (int ig = 0; ig < pw_rho_->npw; ++ig) {
        const int isz = pw_rho_->ig2isz[ig];
        const int iz = isz % pw_rho_->nz;
        const int is = isz / pw_rho_->nz;
        const int ixy = pw_rho_->is2fftixy[is];
        const int ix = ixy / pw_rho_->fftny;
        const int iy = ixy % pw_rho_->fftny;
        ig_of_cell[(ix * pw_rho_->ny + iy) * pw_rho_->nz + iz] = ig;
    }

    std::vector<std::complex<double>> a_r(pw_rho_->nrxx, std::complex<double>(0.0, 0.0));
    std::vector<std::complex<double>> u_r(pw_rho_->nrxx);
    std::vector<std::complex<double>> d_r(pw_rho_->nrxx);
    std::vector<std::complex<double>> d_recip(pw_rho_->npw, std::complex<double>(0.0, 0.0));
    DFPT_KQ_Basis kq;
    for (int ik = 0; ik < nk; ++ik) {
        kq.init(pw_wfc_, q_cart, ik);
        const int npw_kq = kq.get_npwk();
        // k+q stick index -> rho-grid ig through the shared FFT cell
        std::vector<int> kq2rho(npw_kq, -1);
        for (int igl = 0; igl < npw_kq; ++igl) {
            const int isz = kq.get_ig2isz(igl);
            const int iz = isz % pw_wfc_->nz;
            const int is = isz / pw_wfc_->nz;
            const int ixy = pw_wfc_->is2fftixy[is];
            const int ix = ixy / pw_wfc_->fftny;
            const int iy = ixy % pw_wfc_->fftny;
            kq2rho[igl] = ig_of_cell[(ix * pw_rho_->ny + iy) * pw_rho_->nz + iz];
        }
        for (int ib = 0; ib < nbands; ++ib) {
            const double w = wg(ik, ib);
            if (!dfpt_band_occupied(wg, ik, ib)) {
                continue; // unoccupied band: no contribution to the density
            }
            // periodic part u_nk(r) on the shared grid (phase-free FFT)
            pw_wfc_->recip2real(&psi(ik, ib, 0), u_r.data(), ik);
            // periodic part du_nk(r): scatter the k+q coefficients onto the
            // rho grid and transform (same convention, so the product is
            // consistent with the u transform)
            std::fill(d_recip.begin(), d_recip.end(), std::complex<double>(0.0, 0.0));
            const std::vector<std::complex<double>> dpsi = data.get_dpsi(q_idx, ik, ib);
            const int nd = std::min(npw_kq, static_cast<int>(dpsi.size()));
            for (int igl = 0; igl < nd; ++igl) {
                if (kq2rho[igl] >= 0) {
                    d_recip[kq2rho[igl]] = dpsi[igl];
                }
            }
            pw_rho_->recip2real(d_recip.data(), d_r.data());
            // same normalization as the GS density accumulation
            // (elecstate_pw.cpp rhoBandK: w1 = wg / omega)
            const double w1 = w / pw_rho_->omega;
            for (int ir = 0; ir < pw_rho_->nrxx; ++ir) {
                a_r[ir] += w1 * std::conj(u_r[ir]) * d_r[ir];
            }
        }
    }

    // Hermitian completion at q = 0: the band loop above stores only the
    // u_n^* du_n piece; the physical (real) response density also contains
    // the du_n u_n^* piece, whose coefficients are conj(a_{-G}). At q = 0
    // both harmonics coincide: the response is Delta rho = 2 Re a(r), so
    // symmetrize the real-space amplitude before the FFT. The resulting
    // coefficients are exactly Hermitian on the sphere, including
    // one-sided sticks whose -G falls outside it. Away from q = 0 the +q
    // harmonic of the response is exactly the one-sided object and no
    // completion applies.
    const bool q_is_zero = (std::abs(q_frac.x) < 1.0e-10
                            && std::abs(q_frac.y) < 1.0e-10
                            && std::abs(q_frac.z) < 1.0e-10);
    if (q_is_zero) {
        for (int ir = 0; ir < pw_rho_->nrxx; ++ir) {
            a_r[ir] = std::complex<double>(2.0 * a_r[ir].real(), 0.0);
        }
    }

    // q-shifted coefficients A_Delta on the rho grid
    std::vector<std::complex<double>> drho_g(pw_rho_->npw);
    pw_rho_->real2recip(a_r.data(), drho_g.data());

    // charge conservation: the Delta = -q harmonic (G+q = 0 component of the
    // response density) must vanish whenever -q falls on a reciprocal
    // lattice vector; for a generic q inside the cell this never triggers
    {
        const ModuleBase::Vector3<double> mq_cart(-q_cart.x, -q_cart.y, -q_cart.z);
        const ModuleBase::Vector3<double> mfrac = mq_cart * recip_matrix_.Inverse();
        const double mr[3] = {std::round(mfrac.x), std::round(mfrac.y), std::round(mfrac.z)};
        if (std::abs(mfrac.x - mr[0]) < 1.0e-6 &&
            std::abs(mfrac.y - mr[1]) < 1.0e-6 &&
            std::abs(mfrac.z - mr[2]) < 1.0e-6)
        {
            const int ix = (static_cast<int>(mr[0]) % pw_rho_->nx + pw_rho_->nx) % pw_rho_->nx;
            const int iy = (static_cast<int>(mr[1]) % pw_rho_->ny + pw_rho_->ny) % pw_rho_->ny;
            const int iz = (static_cast<int>(mr[2]) % pw_rho_->nz + pw_rho_->nz) % pw_rho_->nz;
            const int ig0 = ig_of_cell[(ix * pw_rho_->ny + iy) * pw_rho_->nz + iz];
            if (ig0 >= 0) {
                drho_g[ig0] = std::complex<double>(0.0, 0.0);
            }
        }
    }
    data.set_drho_g(q_idx, 0, drho_g);

    // real-space manifest density: at q = 0 the completed coefficients are
    // already the full (real) response; away from q = 0 the manifest is the
    // real combination 2 Re[e^{i q r} A(r)] of the one-sided amplitude
    std::vector<std::complex<double>> a_clean(pw_rho_->nrxx);
    pw_rho_->recip2real(drho_g.data(), a_clean.data());
    std::vector<double> drho_r(pw_rho_->nrxx);
    if (q_is_zero) {
        for (int ir = 0; ir < pw_rho_->nrxx; ++ir) {
            drho_r[ir] = a_clean[ir].real();
        }
    } else {
        for (int ix = 0; ix < pw_rho_->nx; ++ix) {
            for (int iy = 0; iy < pw_rho_->ny; ++iy) {
                for (int iz = 0; iz < pw_rho_->nz; ++iz) {
                    const int ir = (ix * pw_rho_->ny + iy) * pw_rho_->nz + iz;
                    const double theta = ModuleBase::TWO_PI *
                                         (q_frac.x * ix / pw_rho_->nx +
                                          q_frac.y * iy / pw_rho_->ny +
                                          q_frac.z * iz / pw_rho_->nz);
                    drho_r[ir] = 2.0 * (a_clean[ir].real() * std::cos(theta) -
                                        a_clean[ir].imag() * std::sin(theta));
                }
            }
        }
    }
    data.set_drho_r(q_idx, 0, drho_r);

    // remember the freshly computed output for the mixing step
    if (q_idx >= static_cast<int>(drho_out_.size())) {
        drho_out_.resize(q_idx + 1);
    }
    drho_out_[q_idx].assign(1, drho_g);
}

void DFPT_Rho::cal_docc(const psi::Psi<std::complex<double>>& psi,
                        const ModuleBase::matrix& wg, int q_idx,
                        DFPT_PW_Data& data) {
    // Reserved first-order occupation matrix (docc) for DFT+U (U0).
    // The physical cross terms need the beta projectors at both k and k+q
    // (a PW-side adapter of the build_vkb machinery); they land together
    // with the Plus_U production wiring in the C7/U1 window, when dpsi and
    // a usable Plus_U provider coexist. Pure-PW runs keep u_active() false
    // and never reach this accumulation:
    //   cross term: Re(becp(k+q, dpsi) * becp(k, psi))      (response)
    //   frozen term: becp(k, psi) * dbecp_f(k, psi)         (GS k)
    if (!data.with_u()) {
        return;
    }
    (void)psi;
    (void)wg;
    (void)q_idx;
    (void)data;
}

void DFPT_Rho::reset_mixing(int q_idx) {
    if (q_idx < 0) {
        return;
    }
    if (q_idx < static_cast<int>(drho_in_.size())) {
        drho_in_[q_idx].clear();
    }
    if (q_idx < static_cast<int>(residual_.size())) {
        residual_[q_idx] = 0.0;
    }
}

void DFPT_Rho::mix_drho(int q_idx, DFPT_PW_Data& data) {
    if (mixer_ == nullptr || pw_rho_ == nullptr) {
        return;
    }
    const std::vector<std::complex<double>> out = data.get_drho_g(q_idx, 0);
    if (out.empty() || static_cast<int>(out.size()) != pw_rho_->npw) {
        return;
    }
    const int npw = pw_rho_->npw;
    if (q_idx >= static_cast<int>(drho_in_.size())) {
        drho_in_.resize(q_idx + 1);
        residual_.resize(q_idx + 1, 0.0);
    }
    // first iteration starts from a zero input density
    if (drho_in_[q_idx].empty()) {
        drho_in_[q_idx].assign(1, std::vector<std::complex<double>>(npw, std::complex<double>(0.0, 0.0)));
    }
    const std::vector<std::complex<double>>& rin = drho_in_[q_idx][0];
    std::vector<std::complex<double>> mixed(npw);
    mixer_->plain_mix(mixed.data(),
                      rin.data(),
                      out.data(),
                      npw,
                      std::function<void(std::complex<double>*)>());
    // relative residual ||out - in|| / ||out||
    double dn2 = 0.0;
    double o2 = 0.0;
    for (int ig = 0; ig < npw; ++ig) {
        dn2 += std::norm(out[ig] - rin[ig]);
        o2 += std::norm(out[ig]);
    }
    residual_[q_idx] = (o2 > 0.0) ? std::sqrt(dn2 / o2) : 0.0;
    drho_in_[q_idx][0] = mixed;
    data.set_drho_g(q_idx, 0, mixed);

    // rebuild the real-space manifest from the mixed coefficients (q = 0:
    // completed coefficients are the full real response; otherwise the
    // one-sided 2 Re[e^{i q r} A(r)] manifest)
    const ModuleBase::Vector3<double> q_frac = data.get_qvec(q_idx);
    const bool q_is_zero = (std::abs(q_frac.x) < 1.0e-10
                            && std::abs(q_frac.y) < 1.0e-10
                            && std::abs(q_frac.z) < 1.0e-10);
    std::vector<std::complex<double>> a_clean(pw_rho_->nrxx);
    pw_rho_->recip2real(mixed.data(), a_clean.data());
    std::vector<double> drho_r(pw_rho_->nrxx);
    if (q_is_zero) {
        for (int ir = 0; ir < pw_rho_->nrxx; ++ir) {
            drho_r[ir] = a_clean[ir].real();
        }
    } else {
        for (int ix = 0; ix < pw_rho_->nx; ++ix) {
            for (int iy = 0; iy < pw_rho_->ny; ++iy) {
                for (int iz = 0; iz < pw_rho_->nz; ++iz) {
                    const int ir = (ix * pw_rho_->ny + iy) * pw_rho_->nz + iz;
                    const double theta = ModuleBase::TWO_PI *
                                         (q_frac.x * ix / pw_rho_->nx +
                                          q_frac.y * iy / pw_rho_->ny +
                                          q_frac.z * iz / pw_rho_->nz);
                    drho_r[ir] = 2.0 * (a_clean[ir].real() * std::cos(theta) -
                                        a_clean[ir].imag() * std::sin(theta));
                }
            }
        }
    }
    data.set_drho_r(q_idx, 0, drho_r);
}

double DFPT_Rho::get_residual(int q_idx, DFPT_PW_Data& data) const {
    (void)data;
    if (q_idx < 0 || q_idx >= static_cast<int>(residual_.size())) {
        return 0.0;
    }
    return residual_[q_idx];
}

void DFPT_Rho::v_hartree_q(const ModuleBase::Vector3<double>& q_cart,
                           const std::vector<std::complex<double>>& drho_g,
                           std::vector<std::complex<double>>& dv_ha_g) const {
    if (pw_rho_ == nullptr) {
        dv_ha_g.clear();
        return;
    }
    const int npw = pw_rho_->npw;
    if (static_cast<int>(drho_g.size()) != npw) {
        dv_ha_g.clear();
        return;
    }
    dv_ha_g.assign(npw, std::complex<double>(0.0, 0.0));
    for (int ig = 0; ig < npw; ++ig) {
        const ModuleBase::Vector3<double> w = pw_rho_->gcar[ig] + q_cart;
        const double w2_lat0 = w * w; // 1/lat0^2 units, like pw_rho_->gg
        // skip |G+q| = 0 (ig = -q): the q-shifted G=0 harmonic of the
        // Hartree kernel (v_hartree skips ig_gge0 the same way)
        if (w2_lat0 < 1.0e-12) {
            continue;
        }
        const double fac = ModuleBase::e2 * ModuleBase::FOUR_PI
                           / (pw_rho_->tpiba2 * w2_lat0);
        dv_ha_g[ig] = fac * drho_g[ig];
    }
}

} // namespace ModuleDFPT
