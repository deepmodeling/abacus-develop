// ============================================================
// DFPT_PW::Impl::solve_displacement and the shared assemble_v_sc
// implementation, with helper extraction.
//
// The two assembly routines (displacement SCF, efield SCF) both need
// the same screened response potential v^{SCF}(drho), so that piece
// is exposed as a separate Impl member reused by both translation
// units via dfpt_pw_impl.h.
// ============================================================

#include "dfpt_pw_impl.h"

#include "dfpt_pw_data.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"

#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace ModuleDFPT
{

void DFPT_PW::Impl::assemble_v_sc(const ModuleBase::Vector3<double>& q_cart,
                                  const std::vector<std::complex<double>>& drho_in_g,
                                  std::vector<std::complex<double>>& v_sc_r) const
{
    const int nrxx = pw_rho_->nrxx;
    v_sc_r.assign(nrxx, std::complex<double>(0.0, 0.0));
    if (drho_in_g.empty() || static_cast<int>(drho_in_g.size()) != pw_rho_->npw)
    {
        return;
    }
    // Hartree q-shifted response: V_H(q,G) = 4pi / |G+q|^2 * drho(G)
    std::vector<std::complex<double>> dv_ha_g;
    rho_.v_hartree_q(q_cart, drho_in_g, dv_ha_g);
    std::vector<std::complex<double>> vh_r(nrxx);
    pw_rho_->recip2real(dv_ha_g.data(), vh_r.data());
    for (int ir = 0; ir < nrxx; ++ir)
    {
        v_sc_r[ir] = vh_r[ir];
    }
    // XC kernel response: v_xc(r) = f_xc(r,rho0) * drho(r); the
    // XC_First_Order provider returns zero-size output for any
    // unsupported / pure-LDA-kernel path, so the size check is the
    // correct existence guard.
    if (xc_ == nullptr)
    {
        return;
    }
    std::vector<std::complex<double>> a_r(nrxx);
    pw_rho_->recip2real(drho_in_g.data(), a_r.data());
    std::vector<std::complex<double>> b_r;
    xc_->apply(a_r, b_r);
    if (static_cast<int>(b_r.size()) != nrxx)
    {
        return;
    }
    for (int ir = 0; ir < nrxx; ++ir)
    {
        v_sc_r[ir] += b_r[ir];
    }
}

double DFPT_PW::Impl::sternheimer_per_band(int ik,
                                           int ib,
                                           const std::vector<std::vector<std::complex<double>>>& dv_sc,
                                           int nbands,
                                           int lin_max,
                                           double lin_thr)
{
    if (!dfpt_band_occupied(wg_, ik, ib))
    {
        return 0.0;
    }
    std::vector<std::complex<double>> rhs = data_.get_dpsi(last_q_, ik, ib);
    const bool dbg = (getenv("DFPT_DEBUG") != nullptr);
    if (rhs.empty() || static_cast<int>(dv_sc.size()) != nbands || rhs.size() != dv_sc[ib].size())
    {
        if (dbg)
        {
            std::cout << "DBG skip solve ik=" << ik << " ib=" << ib << " rhs.size=" << rhs.size()
                      << " dv_sc.size=" << dv_sc.size() << " dv_sc[ib].size="
                      << (dv_sc.size() > static_cast<size_t>(ib) ? dv_sc[ib].size() : 999999) << std::endl;
        }
        return 0.0;
    }
    for (size_t i = 0; i < rhs.size(); ++i)
    {
        rhs[i] = -(rhs[i] + dv_sc[ib][i]);
    }
    hamilt_->set_shift(eig_(ik, ib));
    std::vector<std::complex<double>> dpsi_out;
    double res = 0.0;
    stern_.solve(*hamilt_, occ_kq_[ik], rhs, lin_max, lin_thr, dpsi_out, res);
    if (dbg)
    {
        double nr = 0.0;
        double nb2 = 0.0;
        for (size_t i = 0; i < dpsi_out.size(); ++i)
        {
            nr += std::norm(dpsi_out[i]);
            nb2 += std::norm(rhs[i]);
        }
        std::cout << "DBG solve ik=" << ik << " ib=" << ib << " eps=" << eig_(ik, ib) << " res=" << res
                  << " |dpsi|=" << std::sqrt(nr) << " |rhs|=" << std::sqrt(nb2)
                  << " finite=" << (std::isfinite(std::sqrt(nr)) ? 1 : 0) << std::endl;
    }
    data_.set_dpsi(last_q_, ik, ib, dpsi_out);
    return res;
}

void DFPT_PW::Impl::stash_converged_disp_response(int q_idx,
                                                  int iat,
                                                  int idir,
                                                  const std::vector<std::complex<double>>& v_sc_r_last,
                                                  int nk,
                                                  int nbands)
{
    data_.set_vsc_r(iat, idir, v_sc_r_last);
    std::vector<std::vector<std::vector<std::complex<double>>>> disp(
        nk,
        std::vector<std::vector<std::complex<double>>>(nbands));
    for (int ik = 0; ik < nk; ++ik)
    {
        for (int ib = 0; ib < nbands; ++ib)
        {
            disp[ik][ib] = data_.get_dpsi(q_idx, ik, ib);
        }
    }
    data_.set_dpsi_disp(iat, idir, disp);
}

namespace
{

void debug_iter_snapshot(int iter,
                         const std::vector<std::complex<double>>& drho_in_g,
                         const std::vector<std::complex<double>>& v_sc_r,
                         int npw,
                         int nrxx)
{
    double dh = 0.0;
    double dv = 0.0;
    for (int ig = 0; ig < npw; ++ig)
    {
        dh += std::norm(drho_in_g[ig]);
    }
    for (int ir = 0; ir < nrxx; ++ir)
    {
        dv += std::norm(v_sc_r[ir]);
    }
    std::cout << "DBG iter=" << iter << " |drho_in_g|=" << std::sqrt(dh) << " |v_sc_r|=" << std::sqrt(dv)
              << std::endl;
}

void debug_h_consistency(DFPT_HamiltShift& h,
                         const std::vector<std::vector<std::complex<double>>>& occ_k,
                         const ModuleBase::matrix& eig,
                         int ikq)
{
    std::cout << "DBG occ_kq nstates=" << occ_k.size() << std::endl;
    for (size_t m = 0; m < occ_k.size(); ++m)
    {
        h.set_shift(0.0);
        std::vector<std::complex<double>> hp(occ_k[m].size());
        h.apply(occ_k[m].data(), hp.data());
        std::complex<double> dot(0.0, 0.0);
        for (size_t i = 0; i < hp.size(); ++i)
        {
            dot += std::conj(occ_k[m][i]) * hp[i];
        }
        std::cout << "DBG   <psi_" << m << "|H|psi_" << m << "> = " << dot.real() << " + i " << dot.imag()
                  << "  (GS eig " << eig(ikq, static_cast<int>(m)) << ")" << std::endl;
        std::cout << "DBG   <psi_" << m << "|T+Vnl|psi_" << m << "> = " << h.debug_t_vnl(occ_k[m]) << std::endl;
        std::cout << "DBG   <psi_" << m << "|V(wfc-path)|psi_" << m << "> = " << h.debug_v_wfc(occ_k[m]) << std::endl;
    }
}

} // namespace

double DFPT_PW::Impl::solve_displacement(int q_idx, int iat, int idir)
{
    ModuleBase::TITLE("DFPT_PW", "solve_displacement");
    ModuleBase::timer::start("DFPT_PW", "solve_displacement");
    if (!wired() || hamilt_ == nullptr)
    {
        ModuleBase::timer::end("DFPT_PW", "solve_displacement");
        return 0.0;
    }
    const ModuleBase::Vector3<double> q_frac = data_.get_qvec(q_idx);
    const ModuleBase::Vector3<double> q_cart = q_frac * ucell_->G;
    const int nrxx = pw_rho_->nrxx;
    const int npw = pw_rho_->npw;
    const int nk = gs_psi_.get_nk();
    const int nbands = gs_psi_.get_nbands();
    const bool dbg = (getenv("DFPT_DEBUG") != nullptr);

    pert_.build_dv(q_idx, iat, idir, data_);
    rho_.reset_mixing(q_idx);
    data_.set_drho_g(q_idx, 0, std::vector<std::complex<double>>(npw, std::complex<double>(0.0, 0.0)));

    const int lin_max = data_.get_max_iter();
    const double lin_thr = data_.get_conv_thr();

    bool converged = false;
    double residual = 0.0;
    std::vector<std::complex<double>> v_sc_r_last;
    for (int iter = 0; iter < max_iter_ && !converged; ++iter)
    {
        std::vector<std::complex<double>> v_sc_r;
        const std::vector<std::complex<double>> drho_in_g = data_.get_drho_g(q_idx, 0);
        assemble_v_sc(q_cart, drho_in_g, v_sc_r);
        if (dbg)
        {
            debug_iter_snapshot(iter, drho_in_g, v_sc_r, npw, nrxx);
        }
        v_sc_r_last = v_sc_r;

        for (int ik = 0; ik < nk; ++ik)
        {
            if (static_cast<int>(occ_kq_.size()) <= ik || occ_kq_[ik].empty())
            {
                if (dbg)
                {
                    std::cout << "DBG skip ik=" << ik << " no occ_kq" << std::endl;
                }
                continue;
            }
            pert_.apply_dv(q_idx, ik, gs_psi_, data_);
            std::vector<std::vector<std::complex<double>>> dv_sc;
            pert_.apply_vr(q_idx, ik, v_sc_r, gs_psi_, q_cart, dv_sc);
            if (ik != last_ik_ || last_q_ != q_idx)
            {
                hamilt_->set_context(q_cart, ik);
                last_ik_ = ik;
                last_q_ = q_idx;
                if (dbg)
                {
                    debug_h_consistency(*hamilt_, occ_kq_[ik], eig_, ikq_of_k_[ik]);
                }
            }
            for (int ib = 0; ib < nbands; ++ib)
            {
                sternheimer_per_band(ik, ib, dv_sc, nbands, lin_max, lin_thr);
            }
        }

        rho_.compute_drho(gs_psi_, wg_, q_idx, data_);
        rho_.mix_drho(q_idx, data_);
        residual = rho_.get_residual(q_idx, data_);
        if (dbg)
        {
            std::cout << "DBG iter=" << iter << " residual=" << residual << " conv_thr=" << conv_thr_ << std::endl;
        }
        converged = (residual < conv_thr_);
    }
    stash_converged_disp_response(q_idx, iat, idir, v_sc_r_last, nk, nbands);

    ModuleBase::timer::end("DFPT_PW", "solve_displacement");
    return residual;
}

} // namespace ModuleDFPT
