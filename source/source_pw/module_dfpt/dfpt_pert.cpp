// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in the design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================

#include "dfpt_pert.h"

#include "source_base/constants.h"
#include "source_base/global_function.h"
#include "source_base/math_integral.h"
#include "source_base/math_sphbes.h"
#include "source_base/truncated_func.h"
#include "source_cell/atom_pseudo.h"
#include "source_cell/atom_spec.h"
#include "source_pw/module_pwdft/stru_fac.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>

namespace ModuleDFPT {

DFPT_Pert::DFPT_Pert() {}

DFPT_Pert::~DFPT_Pert() {}

void DFPT_Pert::init(UnitCell& ucell, ModulePW::PW_Basis* pw_rho, 
                     ModulePW::PW_Basis_K* pw_wfc, Structure_Factor& sf) {
    ucell_ = &ucell;
    pw_rho_ = pw_rho;
    pw_wfc_ = pw_wfc;
    sf_ = &sf;
}

void DFPT_Pert::atom_index(int atom_idx, int& it, int& ia) const {
    it = 0;
    ia = atom_idx;
    for (int it_type = 0; it_type < ucell_->ntype; ++it_type) {
        if (ia < ucell_->atoms[it_type].na) {
            it = it_type;
            return;
        }
        ia -= ucell_->atoms[it_type].na;
    }
    // out of range: leave it/ia at the last type / last picture and let the
    // caller guard; dV requests with invalid indices simply produce nothing.
    ia = -1;
}

void DFPT_Pert::rho_gvec(int ig, ModuleBase::Vector3<double>& gcar) const {
    const int isz = pw_rho_->ig2isz[ig];
    int iz = isz % pw_rho_->nz;
    const int is = isz / pw_rho_->nz;
    const int ixy = pw_rho_->is2fftixy[is];
    int ix = ixy / pw_rho_->fftny;
    int iy = ixy % pw_rho_->fftny;
    if (ix >= int(pw_rho_->nx / 2) + 1) { ix -= pw_rho_->nx; }
    if (iy >= int(pw_rho_->ny / 2) + 1) { iy -= pw_rho_->ny; }
    if (iz >= int(pw_rho_->nz / 2) + 1) { iz -= pw_rho_->nz; }
    gcar = ModuleBase::Vector3<double>(ix, iy, iz) * ucell_->G;
}

double DFPT_Pert::vloc_at_g(int it, double g2) const {
    // g2 is the squared magnitude in bohr^-2 units.
    const Atom* atom = &ucell_->atoms[it];
    const double zv = atom->ncpp.zv;
    if (atom->coulomb_potential) {
        // analytic Coulomb local potential (vl_pw.cpp::vloc_coulomb)
        return -zv * ModuleBase::e2 * ModuleBase::FOUR_PI / ucell_->omega / g2;
    }
    // numeric pseudopotential: mirror vl_pw.cpp::vloc_of_g at the requested
    // magnitude instead of interpolating the precomputed shell table. This
    // keeps the rho-grid kernel consistent with the ground-state local
    // potential for every magnitude |Delta+q|.
    const int msh = atom->ncpp.msh;
    const double fac = zv * ModuleBase::e2;
    std::vector<double> aux(msh);
    const double g = std::sqrt(g2);
    if (g < 1.0e-8) {
        double v0 = 0.0;
        for (int ir = 0; ir < msh; ++ir) {
            aux[ir] = atom->ncpp.r[ir] * (atom->ncpp.r[ir] * atom->ncpp.vloc_at[ir] + fac);
        }
        ModuleBase::Integral::Simpson_Integral(msh, aux.data(), atom->ncpp.rab.data(), v0);
        return v0 * ModuleBase::FOUR_PI / ucell_->omega;
    }
    for (int ir = 0; ir < msh; ++ir) {
        aux[ir] = (atom->ncpp.r[ir] * atom->ncpp.vloc_at[ir] + fac * std::erf(atom->ncpp.r[ir]))
                  * std::sin(g * atom->ncpp.r[ir]) / g;
    }
    double v = 0.0;
    ModuleBase::Integral::Simpson_Integral(msh, aux.data(), atom->ncpp.rab.data(), v);
    // erf(r)-compensating gaussian subtraction (same as vloc_of_g)
    v -= fac * ModuleBase::truncated_exp(-g2 * 0.25) / g2;
    return v * ModuleBase::FOUR_PI / ucell_->omega;
}

void DFPT_Pert::dVloc_dtau(int atom_idx, int dir,
                           const ModuleBase::Vector3<double>& q,
                           std::vector<std::complex<double>>& dv) {
    if (pw_rho_ == nullptr || pw_rho_->gamma_only) {
        ModuleBase::WARNING_QUIT("DFPT_Pert::dVloc_dtau",
                                 "DFPT requires a complex (gamma_only=false) real-space basis.");
    }
    int it = 0;
    int ia = 0;
    atom_index(atom_idx, it, ia);
    if (ia < 0) {
        return;
    }
    const ModuleBase::Vector3<double>& tau = ucell_->atoms[it].tau[ia];
    const int npw = pw_rho_->npw;
    dv.assign(npw, std::complex<double>(0.0, 0.0));
    ModuleBase::Vector3<double> gcar;
    for (int ig = 0; ig < npw; ++ig) {
        rho_gvec(ig, gcar);
        const ModuleBase::Vector3<double> w = gcar + q; // Delta + q, 2*pi/lat0 units
        const double w2 = w * w;
        // the Delta == -q component carries no displacement gradient (constant
        // potential shift) and is dropped, consistently with the G=0 handling
        // of the ground-state local potential.
        if (w2 < 1.0e-12) {
            continue;
        }
        const double g_bohr2 = w2 * ucell_->tpiba2;
        const double vloc = vloc_at_g(it, g_bohr2);
        // GS structure-factor convention (stru_fac.cpp: ci_tpi =
        // NEG_IMAG_UNIT * 2pi): exp(-i 2pi (g.tau)), with g in 1/lat0 units
        // and tau in lat0 units; 2pi/lat0 = tpiba only multiplies the
        // magnitude (vl_pw.cpp: qnorm = |g| * tpiba).
        const double arg = -ModuleBase::TWO_PI * (w * tau);
        const std::complex<double> phase(std::cos(arg), std::sin(arg));
        // dV_loc / d tau_direction = -i g_dir * Vloc * exp(-i (Delta+q).tau)
        const std::complex<double> iw_dir =
            std::complex<double>(0.0, -1.0) * (ucell_->tpiba * w[dir]);
        dv[ig] = iw_dir * vloc * phase;
    }
}

void DFPT_Pert::build_dv(int q_idx, int atom_idx, int dir, DFPT_PW_Data& data) {
    // the local first-order potential is assembled on the rho grid in reciprocal
    // space (reciprocal coefficients indexed by the rho-basis ig), then brought
    // to the shared real-space grid where apply_dv performs the convolution.
    if (pw_rho_ == nullptr) {
        return;
    }
    const ModuleBase::Vector3<double> q_cart = data.get_qvec(q_idx) * ucell_->G;
    std::vector<std::complex<double>> dv_recip;
    dVloc_dtau(atom_idx, dir, q_cart, dv_recip);
    data.set_dv_recip_c(q_idx, 0, dv_recip);

    std::vector<std::complex<double>> dv_real(pw_rho_->nrxx);
    pw_rho_->recip2real(dv_recip.data(), dv_real.data());
    data.set_dv_rc(q_idx, 0, dv_real);
    data.set_pert_atom(atom_idx);
    data.set_pert_dir(dir);

    // DFT+U perturbation reservation (U0): append the first-order Hubbard
    // potential dV_U when a DFT+U provider is wired. Physical implementation
    // lands in C1 (frozen projector term) and C3 (occupation response).
    if (data.with_u()) {
        build_dv_u(q_idx, atom_idx, dir, data);
    }
}

void DFPT_Pert::real_space_dv(int q_idx, int k_idx,
                              const psi::Psi<std::complex<double>>& psi,
                              DFPT_PW_Data& data,
                              const DFPT_KQ_Basis& kq,
                              std::vector<std::vector<std::complex<double>>>& dv_psi) const {
    const std::vector<std::complex<double>> dv_rc = data.get_dv_rc(q_idx, 0);
    if (dv_rc.empty() || dv_rc.size() != static_cast<size_t>(pw_rho_->nrxx)) {
        return;
    }
    apply_vr_core(k_idx, dv_rc, psi, kq, dv_psi);
}

void DFPT_Pert::apply_vr(int q_idx, int k_idx,
                         const std::vector<std::complex<double>>& v_rc,
                         const psi::Psi<std::complex<double>>& psi,
                         const ModuleBase::Vector3<double>& q_cart,
                         std::vector<std::vector<std::complex<double>>>& dv_psi) const {
    (void)q_idx;
    if (pw_rho_ == nullptr || pw_wfc_ == nullptr
        || v_rc.size() != static_cast<size_t>(pw_rho_->nrxx)) {
        dv_psi.clear();
        return;
    }
    DFPT_KQ_Basis kq;
    kq.init(pw_wfc_, pw_rho_, q_cart, k_idx);
    apply_vr_core(k_idx, v_rc, psi, kq, dv_psi);
}

void DFPT_Pert::apply_vr_core(int k_idx,
                              const std::vector<std::complex<double>>& v_rc,
                              const psi::Psi<std::complex<double>>& psi,
                              const DFPT_KQ_Basis& kq,
                              std::vector<std::vector<std::complex<double>>>& dv_psi) const {
    const int nbands = psi.get_nbands();
    const int npwk_kq = kq.get_npwk();
    std::vector<std::complex<double>> u_r(pw_rho_->nrxx);
    std::vector<std::complex<double>> d_r(pw_rho_->nrxx);
    std::vector<std::complex<double>> d_recip(pw_rho_->npw);
    dv_psi.assign(nbands, std::vector<std::complex<double>>(npwk_kq, std::complex<double>(0.0, 0.0)));
    for (int iband = 0; iband < nbands; ++iband) {
        pw_wfc_->recip2real(&psi(k_idx, iband, 0), u_r.data(), k_idx);
        for (int ir = 0; ir < pw_rho_->nrxx; ++ir) {
            d_r[ir] = u_r[ir] * v_rc[ir];
        }
        pw_rho_->real2recip(d_r.data(), d_recip.data());
        std::vector<std::complex<double>> dpsi(npwk_kq, std::complex<double>(0.0, 0.0));
        for (int igl = 0; igl < npwk_kq; ++igl) {
            const int ig_rho = kq.get_ig_rho(igl);
            if (ig_rho >= 0) {
                dpsi[igl] = d_recip[ig_rho];
            }
        }
         dv_psi[iband] = dpsi;
     }
}

void DFPT_Pert::apply_dv(int q_idx, int k_idx, const psi::Psi<std::complex<double>>& psi, 
                         DFPT_PW_Data& data) {
    const int atom_idx = data.get_pert_atom();
    const int dir = data.get_pert_dir();
    const ModuleBase::Vector3<double> q_cart = data.get_qvec(q_idx) * ucell_->G;

    DFPT_KQ_Basis kq;
    kq.init(pw_wfc_, pw_rho_, q_cart, k_idx);

    const int nbands = psi.get_nbands();
    std::vector<std::vector<std::complex<double>>> dv_psi(nbands);

    // local contribution: dVloc(r) * psi on the shared real-space grid
    real_space_dv(q_idx, k_idx, psi, data, kq, dv_psi);

    // nonlocal contribution: dVnl/dtau |psi> (per displaced atom)
    std::vector<std::vector<std::complex<double>>> dv_psi_nl;
    dVnl_dtau(atom_idx, dir, q_cart, psi, k_idx, dv_psi_nl);
    if (dv_psi_nl.size() == static_cast<size_t>(nbands)) {
        for (int iband = 0; iband < nbands; ++iband) {
            if (dv_psi[iband].size() != dv_psi_nl[iband].size()) {
                continue;
            }
            for (size_t i = 0; i < dv_psi[iband].size(); ++i) {
                dv_psi[iband][i] += dv_psi_nl[iband][i];
            }
        }
    }

    for (int iband = 0; iband < nbands; ++iband) {
        data.set_dpsi(q_idx, k_idx, iband, dv_psi[iband]);
    }
}

// ---------------------------------------------------------------------------
// nonlocal first-order potential (normal-conserving separable case)
// ---------------------------------------------------------------------------

double DFPT_Pert::radial_vq(int it, int ib, double g) const {
    const pseudo& ncpp = ucell_->atoms[it].ncpp;
    const int l = ncpp.lll[ib];
    int kkbeta = ncpp.kkbeta;
    if (kkbeta > 0 && (kkbeta % 2 == 0)) {
        --kkbeta;
    }
    std::vector<double> jl(kkbeta);
    std::vector<double> aux(kkbeta);
    ModuleBase::Sphbes::Spherical_Bessel(kkbeta, ncpp.r.data(), g, l, jl.data());
    for (int ir = 0; ir < kkbeta; ++ir) {
        aux[ir] = ncpp.betar(ib, ir) * jl[ir] * ncpp.r[ir];
    }
    double v = 0.0;
    ModuleBase::Integral::Simpson_Integral(kkbeta, aux.data(), ncpp.rab.data(), v);
    // tab convention from vnl_pw.cpp: (4pi/sqrt(Omega)) * integral
    return v * ModuleBase::FOUR_PI / std::sqrt(ucell_->omega);
}

double DFPT_Pert::real_ylm(int l, int m, const ModuleBase::Vector3<double>& ghat) const {
    // orthonormal real spherical harmonics Y_{l,m} for l <= 2 with the
    // standard convention, m in [-l, l]:
    //   Y_{l,0}      = sqrt((2l+1)/4pi) P_l^0(cos0)
    //   Y_{l,m>0}    = sqrt(2 (2l+1)/4pi (l-m)!/(l+m)!) P_l^m(cos0) cos(m phi)
    //   Y_{l,m<0}    = sqrt(2 (2l+1)/4pi (l-|m|)!/(l+|m|)!) P_l^{|m|}(cos0) sin(|m| phi)
    // with the associated Legendre convention P_1^1 = -sin0, P_2^1 = -3 sin0 cos0,
    // P_2^2 = 3 sin^2 0. The ABACUS GS vkb applies an additional (-1)^|m| phase
    // for the m>0 channels; exact GS parity is reconciled in the diamond
    // end-to-end test (C7), while the C1 identity test is convention-independent.
    const double x = ghat.x;
    const double y = ghat.y;
    const double z = ghat.z;
    const double r = std::sqrt(x * x + y * y + z * z);
    if (r < 1.0e-12) {
        return (l == 0) ? 0.5 * std::sqrt(1.0 / ModuleBase::PI) : 0.0;
    }
    const double nx = x / r;
    const double ny = y / r;
    const double nz = z / r;
    switch (l) {
        case 0: {
            return 0.5 * std::sqrt(1.0 / ModuleBase::PI);
        }
        case 1: {
            switch (m) {
                case -1: return -0.5 * std::sqrt(3.0 / ModuleBase::PI) * ny;
                case 0: return 0.5 * std::sqrt(3.0 / ModuleBase::PI) * nz;
                case 1: return -0.5 * std::sqrt(3.0 / ModuleBase::PI) * nx;
            }
            break;
        }
        case 2: {
            switch (m) {
                case -2: return 0.5 * std::sqrt(15.0 / ModuleBase::PI) * nx * ny;
                case -1: return -0.5 * std::sqrt(15.0 / ModuleBase::PI) * nz * ny;
                case 0: return 0.25 * std::sqrt(5.0 / ModuleBase::PI) * (3.0 * nz * nz - 1.0);
                case 1: return -0.5 * std::sqrt(15.0 / ModuleBase::PI) * nz * nx;
                case 2: return 0.25 * std::sqrt(15.0 / ModuleBase::PI) * (nx * nx - ny * ny);
            }
            break;
        }
        default: {
            ModuleBase::WARNING_QUIT("DFPT_Pert::real_ylm",
                                     "real_ylm implemented for l<=2 only (DFPT NC path).");
        }
    }
    return 0.0;
}

void DFPT_Pert::build_vkb(int it, int ia,
                          const std::vector<ModuleBase::Vector3<double>>& gk,
                          std::vector<std::vector<std::complex<double>>>& vkb) const {
    // per-type projector bookkeeping mirrors the ground-state vnl_pw.cpp layout:
    // every radial beta (nbeta) with angular momentum l spins out (2l+1)
    // projectors with combined index lm = l^2 + m, m in 0..2l (i.e. the real
    // harmonic m channels -l..l walked as m' = (-1)^(m+1) ceil... ABACUS ylm
    // block: m=0, +1, -1, +2, -2, ...). We use the signed m' directly.
    const pseudo& ncpp = ucell_->atoms[it].ncpp;
    const int nh = ncpp.nh;
    const int ngk = static_cast<int>(gk.size());
    const ModuleBase::Vector3<double>& tau = ucell_->atoms[it].tau[ia];
    vkb.assign(nh, std::vector<std::complex<double>>(ngk, std::complex<double>(0.0, 0.0)));
    if (nh == 0) {
        return;
    }
    int mu = 0;
    for (int ib = 0; ib < ncpp.nbeta; ++ib) {
        const int l = ncpp.lll[ib];
        if (l > 2) {
            ModuleBase::WARNING_QUIT("DFPT_Pert::build_vkb",
                                     "DFPT NC projector path implemented for l<=2 only.");
        }
        const std::complex<double> pref =
            std::pow(std::complex<double>(0.0, -1.0), l); // (-i)^l
        for (int m = 0; m < 2 * l + 1; ++m) {
            // ABACUS real-harmonic walk over the m channels of this radial beta:
            // m=0 -> m'=0; m=1 -> m'=+1; m=2 -> m'=-1; m=3 -> m'=+2; m=4 -> m'=-2
            const int mr = (m == 0) ? 0 : ((m % 2 == 1) ? (m + 1) / 2 : -(m / 2));
            for (int ig = 0; ig < ngk; ++ig) {
                const ModuleBase::Vector3<double>& G = gk[ig]; // k(+q)+G, 2*pi/lat0
                const double gnorm = std::sqrt(G * G) * ucell_->tpiba; // bohr^-1
                // real_ylm handles the |G|=0 point itself (Y_00 is
                // direction-independent; l>0 channels vanish there together
                // with vq), so the raw vector is passed directly.
                const double ylm = real_ylm(l, mr, G);
                const double vq = radial_vq(it, ib, gnorm);
                // GS structure-factor convention (stru_fac.cpp get_sk /
                // eigts, ci_tpi = -2pi i): exp(-i 2pi (gk.tau))
                const double arg = -ModuleBase::TWO_PI * (G * tau);
                const std::complex<double> phase(std::cos(arg), std::sin(arg));
                vkb[mu][ig] = pref * ylm * vq * phase;
            }
            ++mu;
        }
    }
}

void DFPT_Pert::grad_real_ylm(int l, int m, const ModuleBase::Vector3<double>& ghat,
                              double grad[3]) const {
    // analytic gradients of the real_ylm polynomials (l <= 2), consistent
    // with the conventions documented above real_ylm
    const double x = ghat.x;
    const double y = ghat.y;
    const double z = ghat.z;
    const double c1 = 0.5 * std::sqrt(3.0 / ModuleBase::PI);
    const double c2 = 0.5 * std::sqrt(15.0 / ModuleBase::PI);
    const double c20 = 0.25 * std::sqrt(5.0 / ModuleBase::PI);
    grad[0] = grad[1] = grad[2] = 0.0;
    switch (l) {
        case 0:
            return;
        case 1:
            switch (m) {
                case -1: grad[1] = -c1; return;
                case 0: grad[2] = c1; return;
                case 1: grad[0] = -c1; return;
            }
            break;
        case 2:
            switch (m) {
                case -2: grad[0] = c2 * y; grad[1] = c2 * x; return;
                case -1: grad[1] = -c2 * z; grad[2] = -c2 * y; return;
                case 0: grad[2] = 6.0 * c20 * z; return;
                case 1: grad[0] = -c2 * z; grad[2] = -c2 * x; return;
                case 2: grad[0] = 2.0 * c20 * x; grad[1] = -2.0 * c20 * y; return;
            }
            break;
        default:
            ModuleBase::WARNING_QUIT("DFPT_Pert::grad_real_ylm",
                                     "grad_real_ylm implemented for l<=2 only (DFPT NC path).");
    }
}

void DFPT_Pert::build_vkb_dk(int it, int ia, int dir,
                             const std::vector<ModuleBase::Vector3<double>>& gk,
                             std::vector<std::vector<std::complex<double>>>& vkb,
                             std::vector<std::vector<std::complex<double>>>& dvkb) const {
    const pseudo& ncpp = ucell_->atoms[it].ncpp;
    const int nh = ncpp.nh;
    const int ngk = static_cast<int>(gk.size());
    const ModuleBase::Vector3<double>& tau = ucell_->atoms[it].tau[ia];
    if (static_cast<int>(vkb.size()) != nh
        || static_cast<int>(vkb[0].size()) != ngk) {
        ModuleBase::WARNING_QUIT("DFPT_Pert::build_vkb_dk",
                                 "vkb must be built on the same gk list first.");
    }
    dvkb.assign(nh, std::vector<std::complex<double>>(ngk, std::complex<double>(0.0, 0.0)));
    if (nh == 0) {
        return;
    }
    const double dg = 1.0e-4; // bohr^-1, radial central-difference step
    int mu = 0;
    for (int ib = 0; ib < ncpp.nbeta; ++ib) {
        const int l = ncpp.lll[ib];
        const std::complex<double> pref =
            std::pow(std::complex<double>(0.0, -1.0), l); // (-i)^l
        for (int m = 0; m < 2 * l + 1; ++m) {
            const int mr = (m == 0) ? 0 : ((m % 2 == 1) ? (m + 1) / 2 : -(m / 2));
            for (int ig = 0; ig < ngk; ++ig) {
                const ModuleBase::Vector3<double>& G = gk[ig];
                const double gmag = std::sqrt(G * G); // 2*pi/lat0 units
                const double gnorm = gmag * ucell_->tpiba; // bohr^-1
                const double vq0 = radial_vq(it, ib, gnorm);
                const double dvq = (radial_vq(it, ib, gnorm + dg)
                                    - radial_vq(it, ib, std::max(0.0, gnorm - dg)))
                                   / (dg * (gnorm > dg ? 2.0 : 1.0));
                const double arg = -ModuleBase::TWO_PI * (G * tau);
                const std::complex<double> phase(std::cos(arg), std::sin(arg));
                const std::complex<double> dphase =
                    std::complex<double>(0.0, -ModuleBase::TWO_PI * tau[dir]) * phase;
                double dy[3] = {0.0, 0.0, 0.0};
                double ylm = 0.0;
                if (gmag > 1.0e-10) {
                    const ModuleBase::Vector3<double> ghat = G * (1.0 / gmag);
                    ylm = real_ylm(l, mr, ghat);
                    grad_real_ylm(l, mr, ghat, dy);
                    const double gdir[3] = {ghat.x, ghat.y, ghat.z};
                    // chain rule dghat/dk_dir = (e_dir - ghat*ghat_dir)/|G|
                    double dylm_dir = 0.0;
                    for (int c = 0; c < 3; ++c) {
                        dylm_dir += dy[c] * ((c == dir ? 1.0 : 0.0) - gdir[c] * gdir[dir]);
                    }
                    dylm_dir /= gmag;
                    // radial chain: dg/dk_dir = tpiba * ghat_dir
                    const double dradial = dvq * ucell_->tpiba * gdir[dir];
                    dvkb[mu][ig] = pref * phase * (dylm_dir * vq0 + ylm * dradial)
                                   + pref * ylm * vq0 * dphase;
                } else {
                    // degenerate |G| = 0: only the l = 0 channel survives
                    // (real_ylm convention); keep only the phase term
                    ylm = (l == 0) ? 0.5 * std::sqrt(1.0 / ModuleBase::PI) : 0.0;
                    dvkb[mu][ig] = pref * ylm * vq0 * dphase;
                }
            }
            ++mu;
        }
    }
}

void DFPT_Pert::dVnl_dtau(int atom_idx, int dir,
                          const ModuleBase::Vector3<double>& q_cart,
                          const psi::Psi<std::complex<double>>& psi, int k_idx,
                          std::vector<std::vector<std::complex<double>>>& dv_psi) {
    int it = 0;
    int ia = 0;
    atom_index(atom_idx, it, ia);
    if (ia < 0) {
        return;
    }
    const pseudo& ncpp = ucell_->atoms[it].ncpp;
    if (ncpp.tvanp || ncpp.has_so) {
        // the separable NC path documented in C1; ultrasoft and spin-orbit
        // projectors are deferred (their D and augmentation have |k+q| shifts
        // that need the USPP machinery).
        ModuleBase::WARNING_QUIT("DFPT_Pert::dVnl_dtau",
                                 "DFPT nonlocal first-order potential is implemented "
                                 "for normal-conserving separable pseudopotentials only.");
    }
    const int nh = ncpp.nh;

    // projector -> (radial beta index, m channel) table, matching build_vkb.
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

    // incoming k basis: G = k + G' (pw_wfc k-basis index)
    const int npwk = pw_wfc_->npwk[k_idx];
    std::vector<ModuleBase::Vector3<double>> gk_in(npwk);
    for (int ig = 0; ig < npwk; ++ig) {
        gk_in[ig] = pw_wfc_->getgpluskcar(k_idx, ig);
    }
    std::vector<std::vector<std::complex<double>>> vkb_in;
    build_vkb(it, ia, gk_in, vkb_in);

    // outgoing k+q basis
    DFPT_KQ_Basis kq;
    kq.init(pw_wfc_, pw_rho_, q_cart, k_idx);
    const int npwk_kq = kq.get_npwk();
    std::vector<ModuleBase::Vector3<double>> gk_out(npwk_kq);
    for (int igl = 0; igl < npwk_kq; ++igl) {
        gk_out[igl] = kq.get_gpluskq(igl);
    }
    std::vector<std::vector<std::complex<double>>> vkb_out;
    build_vkb(it, ia, gk_out, vkb_out);

    const int nbands = psi.get_nbands();
    dv_psi.assign(nbands, std::vector<std::complex<double>>(npwk_kq, std::complex<double>(0.0, 0.0)));

    for (int iband = 0; iband < nbands; ++iband) {
        // becp_nu(k) = sum_G' conj(vkb_in[nu][G']) psi(G')
        std::vector<std::complex<double>> becp(nh, std::complex<double>(0.0, 0.0));
        for (int nu = 0; nu < nh; ++nu) {
            for (int ig = 0; ig < npwk; ++ig) {
                becp[nu] += std::conj(vkb_in[nu][ig]) * psi(k_idx, iband, ig);
            }
        }
        // dcbecp = D * becp with D_{mu,nu} = dion(ib_mu, ib_nu) delta_{m_mu, m_nu}
        std::vector<std::complex<double>> dcbecp(nh, std::complex<double>(0.0, 0.0));
        for (int mu = 0; mu < nh; ++mu) {
            for (int nu = 0; nu < nh; ++nu) {
                if (mu_m[mu] != mu_m[nu]) {
                    continue;
                }
                dcbecp[mu] += ncpp.dion(mu_ib[mu], mu_ib[nu]) * becp[nu];
            }
        }
        // term A: i (k+q+G'')_dir * (Vnl |psi>) on the k+q basis
        std::vector<std::complex<double>> term_a(npwk_kq, std::complex<double>(0.0, 0.0));
        for (int igl = 0; igl < npwk_kq; ++igl) {
            std::complex<double> vnlpsi(0.0, 0.0);
            for (int mu = 0; mu < nh; ++mu) {
                vnlpsi += vkb_out[mu][igl] * dcbecp[mu];
            }
            term_a[igl] = std::complex<double>(0.0, 1.0) * (ucell_->tpiba * gk_out[igl][dir]) * vnlpsi;
        }
        // term B: Vnl [i (k+G')_dir |psi>]
        std::vector<std::complex<double>> becp_dpsi(nh, std::complex<double>(0.0, 0.0));
        for (int nu = 0; nu < nh; ++nu) {
            for (int ig = 0; ig < npwk; ++ig) {
                const std::complex<double> dpsi_ig =
                    std::complex<double>(0.0, 1.0) * (ucell_->tpiba * gk_in[ig][dir]) * psi(k_idx, iband, ig);
                becp_dpsi[nu] += std::conj(vkb_in[nu][ig]) * dpsi_ig;
            }
        }
        std::vector<std::complex<double>> dcbecp_dpsi(nh, std::complex<double>(0.0, 0.0));
        for (int mu = 0; mu < nh; ++mu) {
            for (int nu = 0; nu < nh; ++nu) {
                if (mu_m[mu] != mu_m[nu]) {
                    continue;
                }
                dcbecp_dpsi[mu] += ncpp.dion(mu_ib[mu], mu_ib[nu]) * becp_dpsi[nu];
            }
        }
        std::vector<std::complex<double>> term_b(npwk_kq, std::complex<double>(0.0, 0.0));
        for (int igl = 0; igl < npwk_kq; ++igl) {
            std::complex<double> vnl_dpsi(0.0, 0.0);
            for (int mu = 0; mu < nh; ++mu) {
                vnl_dpsi += vkb_out[mu][igl] * dcbecp_dpsi[mu];
            }
            term_b[igl] = vnl_dpsi;
        }
        for (int igl = 0; igl < npwk_kq; ++igl) {
            // GS exp(-2pi gk.tau) projector convention: dVnl/dtau_dir
            // |psi> = -i (k+q+G'')_dir (Vnl|psi>) + Vnl[i (k+G')_dir |psi>]
            dv_psi[iband][igl] = term_b[igl] - term_a[igl];
        }
    }
}

void DFPT_Pert::build_dv_u(int q_idx, int atom_idx, int dir, DFPT_PW_Data& data) {
    // C1 frozen term of the first-order Hubbard potential:
    //   |dphi(k+q)/dtau> V_eff <phi(k)|psi> + adjoint
    // The provider is only usable when its occupation matrices are
    // initialized (u_active()); DFPT_PW::init additionally rejects a wired
    // provider outright (every U hook is a no-op U0 reservation), so this
    // guard is defense in depth; the diamond DFT+U test (C7) will exercise
    // this path once OnsiteProjector integration on the DFPT k+q basis is
    // finalized.
    if (!data.u_active()) {
        return;
    }
    (void)q_idx;
    (void)atom_idx;
    (void)dir;
    // TODO(C7): Onsite_Proj_tools on the DFPT k+q basis; occupation response
    // (|phi(k+q)> U(diag*delta - docc) <phi(k)|psi>) lands in C3 after docc.
}

void DFPT_Pert::d2vloc_r(int atom_idx, int da, int db,
                         std::vector<std::complex<double>>& dv2_r) const {
    if (pw_rho_ == nullptr) {
        return;
    }
    int it = 0;
    int ia = 0;
    atom_index(atom_idx, it, ia);
    if (ia < 0) {
        dv2_r.clear();
        return;
    }
    const ModuleBase::Vector3<double>& tau = ucell_->atoms[it].tau[ia];
    const int npw = pw_rho_->npw;
    std::vector<std::complex<double>> dv2_recip(npw, std::complex<double>(0.0, 0.0));
    ModuleBase::Vector3<double> gcar;
    for (int ig = 0; ig < npw; ++ig) {
        rho_gvec(ig, gcar);
        // QE ground truth (dynmat_us.f90): the mixed (+q,-q) second-order
        // local potential is the integer-G, q-independent kernel
        // -tpiba^2 G_da G_db vloc(|G|) exp(-i G.tau); the (+q,-q) dressings
        // collapse the carrier to 0 for every q, not only when 2q is
        // reciprocal.
        const ModuleBase::Vector3<double> w = gcar;
        const double w2 = w * w;
        if (w2 < 1.0e-12) {
            continue;
        }
        const double vloc = vloc_at_g(it, w2 * ucell_->tpiba2);
        const double arg = -ModuleBase::TWO_PI * (w * tau);
        const std::complex<double> phase(std::cos(arg), std::sin(arg));
        // (-i g_da)(-i g_db) = -g_da g_db
        dv2_recip[ig] = -(ucell_->tpiba * w[da]) * (ucell_->tpiba * w[db]) * vloc * phase;
    }
    dv2_r.assign(pw_rho_->nrxx, std::complex<double>(0.0, 0.0));
    pw_rho_->recip2real(dv2_recip.data(), dv2_r.data());
}

void DFPT_Pert::apply_d2vnl(int atom_idx, int da, int db,
                            const ModuleBase::Vector3<double>& q_eff,
                            bool include_middle,
                            const psi::Psi<std::complex<double>>& psi, int k_idx,
                            std::vector<std::vector<std::complex<double>>>& d2v_psi) const {
    int it = 0;
    int ia = 0;
    atom_index(atom_idx, it, ia);
    if (ia < 0) {
        return;
    }
    const pseudo& ncpp = ucell_->atoms[it].ncpp;
    if (ncpp.tvanp || ncpp.has_so) {
        ModuleBase::WARNING_QUIT("DFPT_Pert::apply_d2vnl",
                                 "DFPT second-order nonlocal potential is implemented "
                                 "for normal-conserving separable pseudopotentials only.");
    }
    const int nh = ncpp.nh;
    const int nbands = psi.get_nbands();

    // projector -> (radial index, m channel) table, matching build_vkb
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

    // incoming k basis and outgoing k+q basis projectors (same atom)
    const int npwk = pw_wfc_->npwk[k_idx];
    std::vector<ModuleBase::Vector3<double>> gk_in(npwk);
    for (int ig = 0; ig < npwk; ++ig) {
        gk_in[ig] = pw_wfc_->getgpluskcar(k_idx, ig);
    }
    std::vector<std::vector<std::complex<double>>> vkb_in;
    build_vkb(it, ia, gk_in, vkb_in);
    DFPT_KQ_Basis kq;
    kq.init(pw_wfc_, pw_rho_, q_eff, k_idx);
    const int npwk_kq = kq.get_npwk();
    std::vector<ModuleBase::Vector3<double>> gk_out(npwk_kq);
    for (int igl = 0; igl < npwk_kq; ++igl) {
        gk_out[igl] = kq.get_gpluskq(igl);
    }
    std::vector<std::vector<std::complex<double>>> vkb_out;
    build_vkb(it, ia, gk_out, vkb_out);

    d2v_psi.assign(nbands, std::vector<std::complex<double>>(npwk_kq, std::complex<double>(0.0, 0.0)));
    for (int iband = 0; iband < nbands; ++iband) {
        // becp and its (k+G')-weighted variants: becp_x = sum x(G') |beta><psi|
        std::vector<std::complex<double>> becp(nh, std::complex<double>(0.0, 0.0));
        std::vector<std::complex<double>> becp_a(nh, std::complex<double>(0.0, 0.0));
        std::vector<std::complex<double>> becp_b(nh, std::complex<double>(0.0, 0.0));
        std::vector<std::complex<double>> becp_ab(nh, std::complex<double>(0.0, 0.0));
        for (int nu = 0; nu < nh; ++nu) {
            for (int ig = 0; ig < npwk; ++ig) {
                const std::complex<double> vc = std::conj(vkb_in[nu][ig]) * psi(k_idx, iband, ig);
                const double kp_da = ucell_->tpiba * gk_in[ig][da];
                const double kp_db = ucell_->tpiba * gk_in[ig][db];
                becp[nu] += vc;
                becp_a[nu] += kp_da * vc;
                becp_b[nu] += kp_db * vc;
                becp_ab[nu] += kp_da * kp_db * vc;
            }
        }
        // D contraction with the same-m selection rule as dVnl_dtau
        std::vector<std::complex<double>> d0(nh, std::complex<double>(0.0, 0.0));
        std::vector<std::complex<double>> da_(nh, std::complex<double>(0.0, 0.0));
        std::vector<std::complex<double>> db_(nh, std::complex<double>(0.0, 0.0));
        std::vector<std::complex<double>> dab(nh, std::complex<double>(0.0, 0.0));
        for (int mu = 0; mu < nh; ++mu) {
            for (int nu = 0; nu < nh; ++nu) {
                if (mu_m[mu] != mu_m[nu]) {
                    continue;
                }
                const double dij = ncpp.dion(mu_ib[mu], mu_ib[nu]);
                d0[mu] += dij * becp[nu];
                da_[mu] += dij * becp_a[nu];
                db_[mu] += dij * becp_b[nu];
                dab[mu] += dij * becp_ab[nu];
            }
        }
        // chi(G'') = sum_mu vkb_out,mu [ -kq_da kq_db d0 - dab
        //                              + (include_middle ? kq_da db_ + kq_db da_ : 0) ]_mu
        // QE ground truth (dynmat_us.f90 + phq_init.f90): the KB second-order
        // term pairs gammap (integer-G (k+G)_da(k+G)_db derivative of beta)
        // with becp1 = <beta_k|psi_k> and the same-atom alphap_a* alphap_b
        // middle product; everything is built at k with integer-G momentum
        // factors, so the caller passes q_eff = 0 and the kernel is
        // q-independent for every q.
        for (int igl = 0; igl < npwk_kq; ++igl) {
            const double kq_da = ucell_->tpiba * gk_out[igl][da];
            const double kq_db = ucell_->tpiba * gk_out[igl][db];
            std::complex<double> chi(0.0, 0.0);
            for (int mu = 0; mu < nh; ++mu) {
                chi += vkb_out[mu][igl] * (-kq_da * kq_db * d0[mu] - dab[mu]);
                if (include_middle) {
                    chi += vkb_out[mu][igl] * (kq_da * db_[mu] + kq_db * da_[mu]);
                }
            }
            d2v_psi[iband][igl] = chi;
        }
    }
}

void DFPT_Pert::build_efield(const ModuleBase::Vector3<double>& field, DFPT_PW_Data& data) {
    // first-order electric-field potential: delta V(r) = - r . E (q=0 limit,
    // position operator in the periodic cell). Computed directly on the shared
    // real-space grid. Only relevant for the Q0 dielectric response (C6).
    if (pw_rho_ == nullptr) {
        return;
    }
    if (pw_rho_->gamma_only) {
        ModuleBase::WARNING_QUIT("DFPT_Pert::build_efield",
                                 "DFPT requires a complex (gamma_only=false) real-space basis.");
    }
    std::vector<std::complex<double>> dv_real(pw_rho_->nrxx, std::complex<double>(0.0, 0.0));
    const ModuleBase::Matrix3& latvec = ucell_->latvec;
    const double lat0 = ucell_->lat0;
    // shared real-space grid layout (serial pool): ir = (ix*ny + iy)*nz + iz,
    // i.e. z runs fastest (verified against the impulse response of the FFT).
    for (int ir = 0; ir < pw_rho_->nrxx; ++ir) {
        const int iz = ir % pw_rho_->nz;
        const int rem = ir / pw_rho_->nz;
        const int iy = rem % pw_rho_->ny;
        const int ix = rem / pw_rho_->ny;
        const double fx = ix / static_cast<double>(pw_rho_->nx);
        const double fy = iy / static_cast<double>(pw_rho_->ny);
        const double fz = iz / static_cast<double>(pw_rho_->nz);
        ModuleBase::Vector3<double> r;
        r.x = (fx * latvec.e11 + fy * latvec.e12 + fz * latvec.e13) * lat0;
        r.y = (fx * latvec.e21 + fy * latvec.e22 + fz * latvec.e23) * lat0;
        r.z = (fx * latvec.e31 + fy * latvec.e32 + fz * latvec.e33) * lat0;
        dv_real[ir] = -(field * r); // -e r.E (e absorbed in field convention)
    }
    data.set_dv_rc(0, 0, dv_real);
}

} // namespace ModuleDFPT