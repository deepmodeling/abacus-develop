// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================

#include "dfpt_hamilt_shift.h"

#include "dfpt_pert.h"
#include "source_base/constants.h"
#include "source_base/global_function.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_basis/module_pw/pw_basis_k.h"
#include "source_cell/unitcell.h"

#include <cmath>
#include <complex>
#include <vector>

namespace ModuleDFPT {

DFPT_HamiltShift::DFPT_HamiltShift(const UnitCell& ucell,
                                   ModulePW::PW_Basis* pw_rho,
                                   ModulePW::PW_Basis_K* pw_wfc,
                                   const std::vector<double>& veff_r,
                                   const DFPT_Pert* pert)
    : ucell_(&ucell),
      pw_rho_(pw_rho),
      pw_wfc_(pw_wfc),
      pert_(pert),
      veff_r_(veff_r),
      tpiba2_(ucell.tpiba2),
      nrxx_(pw_rho != nullptr ? pw_rho->nrxx : 0) {
    for (int it = 0; it < ucell_->ntype; ++it) {
        const pseudo& ncpp = ucell_->atoms[it].ncpp;
        if (ncpp.tvanp || ncpp.has_so) {
            ModuleBase::WARNING_QUIT("DFPT_HamiltShift",
                                     "the shifted Sternheimer operator is implemented for "
                                     "normal-conserving separable pseudopotentials only.");
        }
        // projector -> (radial beta index, m channel) tables, matching
        // build_vkb / dVnl_dtau
        std::vector<int> ib;
        std::vector<int> m;
        int mu = 0;
        for (int ibeta = 0; ibeta < ncpp.nbeta; ++ibeta) {
            const int l = ncpp.lll[ibeta];
            for (int im = 0; im < 2 * l + 1; ++im) {
                if (mu < ncpp.nh) {
                    ib.push_back(ibeta);
                    m.push_back(im);
                }
                ++mu;
            }
        }
        mu_ib_.push_back(ib);
        mu_m_.push_back(m);
    }
}

DFPT_HamiltShift::~DFPT_HamiltShift() {}

void DFPT_HamiltShift::set_context(const ModuleBase::Vector3<double>& q_cart, int k_idx) {
    kq_.init(pw_wfc_, q_cart, k_idx);
    ik_cache_ = k_idx;
    const int npw = kq_.get_npwk();

    // rho ig -> shared FFT-cell reverse map, then k+q -> rho through the
    // (ix,iy,iz) triple (the stick encodings of the two bases differ)
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
    kq2rho_.assign(npw, -1);
    for (int igl = 0; igl < npw; ++igl) {
        const int isz = kq_.get_ig2isz(igl);
        const int iz = isz % pw_wfc_->nz;
        const int is = isz / pw_wfc_->nz;
        const int ixy = pw_wfc_->is2fftixy[is];
        const int ix = ixy / pw_wfc_->fftny;
        const int iy = ixy % pw_wfc_->fftny;
        kq2rho_[igl] = ig_of_cell[(ix * pw_rho_->ny + iy) * pw_rho_->nz + iz];
    }

    // cache the beta projectors of every atom on the k+q list
    std::vector<ModuleBase::Vector3<double>> gk(npw);
    for (int igl = 0; igl < npw; ++igl) {
        gk[igl] = kq_.get_gpluskq(igl);
    }
    vkb_cache_.assign(ucell_->nat, std::vector<std::vector<std::complex<double>>>());
    for (int iat = 0; iat < ucell_->nat; ++iat) {
        const int it = ucell_->iat2it[iat];
        const int ia = ucell_->iat2ia[iat];
        if (ucell_->atoms[it].ncpp.nh == 0) {
            continue;
        }
        pert_->build_vkb(it, ia, gk, vkb_cache_[iat]);
    }

    x_recip_.assign(pw_rho_->npw, std::complex<double>(0.0, 0.0));
    x_r_.assign(nrxx_, std::complex<double>(0.0, 0.0));
}

void DFPT_HamiltShift::set_shift(double shift) {
    shift_ = shift;
}

int DFPT_HamiltShift::dimension() const {
    return kq_.get_npwk();
}

void DFPT_HamiltShift::apply(const std::complex<double>* x, std::complex<double>* y) const {
    const int npw = kq_.get_npwk();
    if (npw <= 0 || x == nullptr || y == nullptr) {
        return;
    }
    // kinetic part minus the eigenvalue shift
    for (int igl = 0; igl < npw; ++igl) {
        y[igl] = (tpiba2_ * kq_.get_gk2(igl) - shift_) * x[igl];
    }
    // local effective potential: phase-free FFT convolution on the shared
    // grid (the k+q Bloch phases cancel in the product, real_space_dv conv.)
    std::fill(x_recip_.begin(), x_recip_.end(), std::complex<double>(0.0, 0.0));
    for (int igl = 0; igl < npw; ++igl) {
        if (kq2rho_[igl] >= 0) {
            x_recip_[kq2rho_[igl]] = x[igl];
        }
    }
    pw_rho_->recip2real(x_recip_.data(), x_r_.data());
    for (int ir = 0; ir < nrxx_; ++ir) {
        x_r_[ir] *= veff_r_[ir];
    }
    pw_rho_->real2recip(x_r_.data(), x_recip_.data());
    for (int igl = 0; igl < npw; ++igl) {
        if (kq2rho_[igl] >= 0) {
            y[igl] += x_recip_[kq2rho_[igl]];
        }
    }
    // nonlocal part with the cached k+q projectors
    for (int iat = 0; iat < ucell_->nat; ++iat) {
        const int it = ucell_->iat2it[iat];
        const int nh = ucell_->atoms[it].ncpp.nh;
        if (nh == 0) {
            continue;
        }
        const std::vector<std::vector<std::complex<double>>>& vkb = vkb_cache_[iat];
        becp_.assign(nh, std::complex<double>(0.0, 0.0));
        for (int mu = 0; mu < nh; ++mu) {
            for (int igl = 0; igl < npw; ++igl) {
                becp_[mu] += std::conj(vkb[mu][igl]) * x[igl];
            }
        }
        dbecp_.assign(nh, std::complex<double>(0.0, 0.0));
        for (int mu = 0; mu < nh; ++mu) {
            for (int nu = 0; nu < nh; ++nu) {
                if (mu_m_[it][mu] != mu_m_[it][nu]) {
                    continue;
                }
                dbecp_[mu] += ucell_->atoms[it].ncpp.dion(mu_ib_[it][mu], mu_ib_[it][nu]) * becp_[nu];
            }
        }
        for (int mu = 0; mu < nh; ++mu) {
            for (int igl = 0; igl < npw; ++igl) {
                y[igl] += vkb[mu][igl] * dbecp_[mu];
            }
        }
    }
}

double DFPT_HamiltShift::debug_t_vnl(const std::vector<std::complex<double>>& x) const {
    const int npw = kq_.get_npwk();
    double ekin = 0.0;
    for (int igl = 0; igl < npw; ++igl) {
        ekin += tpiba2_ * kq_.get_gk2(igl) * std::norm(x[igl]);
    }
    double vnl = 0.0;
    for (int iat = 0; iat < ucell_->nat; ++iat) {
        const int it = ucell_->iat2it[iat];
        const int nh = ucell_->atoms[it].ncpp.nh;
        if (nh == 0) {
            continue;
        }
        const std::vector<std::vector<std::complex<double>>>& vkb = vkb_cache_[iat];
        becp_.assign(nh, std::complex<double>(0.0, 0.0));
        for (int mu = 0; mu < nh; ++mu) {
            for (int igl = 0; igl < npw; ++igl) {
                becp_[mu] += std::conj(vkb[mu][igl]) * x[igl];
            }
        }
        dbecp_.assign(nh, std::complex<double>(0.0, 0.0));
        for (int mu = 0; mu < nh; ++mu) {
            for (int nu = 0; nu < nh; ++nu) {
                if (mu_m_[it][mu] != mu_m_[it][nu]) {
                    continue;
                }
                dbecp_[mu] += ucell_->atoms[it].ncpp.dion(mu_ib_[it][mu], mu_ib_[it][nu]) * becp_[nu];
            }
        }
        for (int mu = 0; mu < nh; ++mu) {
            vnl += std::real(std::conj(becp_[mu]) * dbecp_[mu]);
        }
    }
    return ekin + vnl;
}

double DFPT_HamiltShift::debug_v_wfc(const std::vector<std::complex<double>>& x) const {
    const int npw = kq_.get_npwk();
    std::vector<std::complex<double>> ur(nrxx_, std::complex<double>(0.0, 0.0));
    pw_wfc_->recip2real(x.data(), ur.data(), ik_cache_);
    for (int ir = 0; ir < nrxx_; ++ir) {
        ur[ir] *= veff_r_[ir];
    }
    std::vector<std::complex<double>> xg(pw_wfc_->npwk[ik_cache_], std::complex<double>(0.0, 0.0));
    pw_wfc_->real2recip(ur.data(), xg.data(), ik_cache_);
    std::complex<double> dot(0.0, 0.0);
    const int n = std::min(static_cast<int>(xg.size()), npw);
    for (int igl = 0; igl < n; ++igl) {
        dot += std::conj(x[igl]) * xg[igl];
    }
    return dot.real();
}

} // namespace ModuleDFPT
