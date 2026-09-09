// ============================================================
// DFPT_PW::Impl::solve_pos_resp and solve_efield_resp
// implementations (q=0 legs). The position response Y^a solves the
// velocity-form commutator Sternheimer equation (Giannozzi et al.
// 1991 QE dvpsi_e), and the e-field SCF response reuses the shared
// Impl::assemble_v_sc from dfpt_pw_solve.cpp. Both large drivers
// are split into small private helpers so the new translation unit
// stays well under the coding-rule 500-line budget and each helper
// function remains under cyclomatic complexity 10.
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

void DFPT_PW::Impl::vel_diag_part(int ik,
                                  int a,
                                  int nbands,
                                  std::vector<std::vector<std::complex<double>>>& vel) const
{
    // dH/dk_a |psi> diagonal kinetic part:
    //   - d/dk_a ( (hbar^2/2m) (k+G)^2 ) = hbar^2 (k+G)_a,
    // re-expressed in ABACUS Rydberg units via tpiba^2.
    const double tpiba2 = ucell_->tpiba * ucell_->tpiba;
    const int npwk = pw_wfc_->npwk[ik];
    vel.assign(nbands, std::vector<std::complex<double>>(npwk, std::complex<double>(0.0, 0.0)));
    for (int ib = 0; ib < nbands; ++ib)
    {
        for (int ig = 0; ig < npwk; ++ig)
        {
            const ModuleBase::Vector3<double> gk = pw_wfc_->getgpluskcar(ik, ig);
            vel[ib][ig] = 2.0 * tpiba2 * gk[a] * gs_psi_(ik, ib, ig);
        }
    }
}

void DFPT_PW::Impl::vel_nl_per_atom(int ik,
                                    int a,
                                    int it,
                                    int ia,
                                    int nbands,
                                    const std::vector<ModuleBase::Vector3<double>>& gk,
                                    std::vector<std::vector<std::complex<double>>>& vel)
{
    const pseudo& ncpp = ucell_->atoms[it].ncpp;
    const int nh = ncpp.nh;
    if (nh == 0)
    {
        return;
    }
    const int npwk = pw_wfc_->npwk[ik];
    // projector -> (radial beta index, m channel) table
    std::vector<int> mu_ib(nh, 0);
    std::vector<int> mu_m(nh, 0);
    int mu_idx = 0;
    for (int ib = 0; ib < ncpp.nbeta; ++ib)
    {
        const int l = ncpp.lll[ib];
        for (int m = 0; m < 2 * l + 1; ++m)
        {
            if (mu_idx < nh)
            {
                mu_ib[mu_idx] = ib;
                mu_m[mu_idx] = m;
            }
            ++mu_idx;
        }
    }
    std::vector<std::vector<std::complex<double>>> vkb;
    pert_.build_vkb(it, ia, gk, vkb);
    // becp_b[mu] = <vkb_mu|psi_b>
    std::vector<std::vector<std::complex<double>>> becp(nbands);
    for (int b = 0; b < nbands; ++b)
    {
        becp[b].assign(nh, std::complex<double>(0.0, 0.0));
        for (int mu = 0; mu < nh; ++mu)
        {
            for (int ig = 0; ig < npwk; ++ig)
            {
                becp[b][mu] += std::conj(vkb[mu][ig]) * gs_psi_(ik, b, ig);
            }
        }
    }
    std::vector<std::vector<std::complex<double>>> dvkb;
    pert_.build_vkb_dk(it, ia, a, gk, vkb, dvkb);
    // dbecp_b[mu] = <dvkb_mu|psi_b>
    std::vector<std::vector<std::complex<double>>> dbecp(nbands);
    for (int b = 0; b < nbands; ++b)
    {
        dbecp[b].assign(nh, std::complex<double>(0.0, 0.0));
        for (int mu = 0; mu < nh; ++mu)
        {
            for (int ig = 0; ig < npwk; ++ig)
            {
                dbecp[b][mu] += std::conj(dvkb[mu][ig]) * gs_psi_(ik, b, ig);
            }
        }
    }
    // dV_nl/dk_a|psi_b> = sum_mu |dvkb_mu> (D becp_b)_mu
    //                              + |vkb_mu> (D dbecp_b)_mu
    for (int b = 0; b < nbands; ++b)
    {
        for (int mu = 0; mu < nh; ++mu)
        {
            std::complex<double> out_b(0.0, 0.0);
            std::complex<double> in_b(0.0, 0.0);
            for (int nu = 0; nu < nh; ++nu)
            {
                if (mu_m[mu] != mu_m[nu])
                {
                    continue;
                }
                const double dij = ncpp.dion(mu_ib[mu], mu_ib[nu]);
                out_b += dij * becp[b][nu];
                in_b += dij * dbecp[b][nu];
            }
            for (int ig = 0; ig < npwk; ++ig)
            {
                vel[b][ig] += dvkb[mu][ig] * out_b + vkb[mu][ig] * in_b;
            }
        }
    }
}

void DFPT_PW::Impl::pos_per_band_solve(int ik,
                                       int a,
                                       int nbands,
                                       int lin_max,
                                       double lin_thr,
                                       std::vector<std::vector<std::vector<std::complex<double>>>>& yvec)
{
    const bool dbg = (getenv("DFPT_DEBUG") != nullptr);
    const double tpiba = ucell_->tpiba;
    const int npwk = pw_wfc_->npwk[ik];
    std::vector<std::vector<std::complex<double>>> vel(nbands);
    vel_diag_part(ik, a, nbands, vel);
    for (int it = 0; it < ucell_->ntype; ++it)
    {
        for (int ia = 0; ia < ucell_->atoms[it].na; ++ia)
        {
            std::vector<ModuleBase::Vector3<double>> gk(npwk);
            for (int ig = 0; ig < npwk; ++ig)
            {
                gk[ig] = pw_wfc_->getgpluskcar(ik, ig);
            }
            vel_nl_per_atom(ik, a, it, ia, nbands, gk, vel);
        }
    }
    // solve (H - eps_v) Y = -(i/tpiba) vel for every occupied band
    for (int ib = 0; ib < nbands; ++ib)
    {
        if (!dfpt_band_occupied(wg_, ik, ib))
        {
            continue;
        }
        std::vector<std::complex<double>> rhs(npwk, std::complex<double>(0.0, 0.0));
        const std::complex<double> fac(0.0, -1.0 / tpiba);
        for (int ig = 0; ig < npwk; ++ig)
        {
            rhs[ig] = fac * vel[ib][ig];
        }
        hamilt_->set_shift(eig_(ik, ib));
        double res = 0.0;
        stern_.solve(*hamilt_, occ_kq_[ik], rhs, lin_max, lin_thr, yvec[ik][ib], res);
        if (dbg)
        {
            std::cout << "DBG posresp a=" << a << " ik=" << ik << " ib=" << ib << " eps=" << eig_(ik, ib)
                      << " res=" << res << std::endl;
        }
    }
}

void DFPT_PW::Impl::efield_per_band_solve(int ik,
                                          int a,
                                          int nbands,
                                          int lin_max,
                                          double lin_thr,
                                          const std::vector<std::vector<std::vector<std::complex<double>>>>& yr,
                                          const std::vector<std::vector<std::complex<double>>>& dv_sc)
{
    (void)a; // used only by caller-index semantics, kept future-proof
    for (int ib = 0; ib < nbands; ++ib)
    {
        if (!dfpt_band_occupied(wg_, ik, ib))
        {
            continue;
        }
        if (yr[ik][ib].empty() || static_cast<int>(dv_sc.size()) != nbands
            || yr[ik][ib].size() != dv_sc[ib].size())
        {
            continue;
        }
        std::vector<std::complex<double>> rhs(yr[ik][ib].size());
        for (size_t i = 0; i < rhs.size(); ++i)
        {
            rhs[i] = -(yr[ik][ib][i] + dv_sc[ib][i]);
        }
        hamilt_->set_shift(eig_(ik, ib));
        std::vector<std::complex<double>> dpsi_out;
        double res = 0.0;
        stern_.solve(*hamilt_, occ_kq_[ik], rhs, lin_max, lin_thr, dpsi_out, res);
        data_.set_dpsi(last_q_, ik, ib, dpsi_out);
        (void)res; // convergence is aggregated via drho residual, not per-band stern_
    }
}

void DFPT_PW::Impl::stash_dpsi_efield(int q_idx, int a, int nk, int nbands)
{
    // stash dpsi^E,a before any later solve reuses the slots
    std::vector<std::vector<std::vector<std::complex<double>>>> de(
        nk,
        std::vector<std::vector<std::complex<double>>>(nbands));
    for (int ik = 0; ik < nk; ++ik)
    {
        for (int ib = 0; ib < nbands; ++ib)
        {
            de[ik][ib] = data_.get_dpsi(q_idx, ik, ib);
        }
    }
    data_.set_dpsi_efield(a, de);
}

void DFPT_PW::Impl::solve_pos_resp(int q_idx)
{
    ModuleBase::TITLE("DFPT_PW", "solve_pos_resp");
    ModuleBase::timer::start("DFPT_PW", "solve_pos_resp");
    // Y^a_{k,v} = P_c x_a|psi_{k,v}> through the Sternheimer equation
    //   (H(k) - eps_v) Y^a_v = P_c [H, x_a]|psi_v>,
    //   [H, x_a]|psi> = -(i/tpiba) dH/dk_a|psi>   (velocity form),
    // exactly the linear solve of QE dvpsi_e (whose rhs negation restores
    // P_c[H,x]psi from commutator_Hx_psi's [x,H] convention). dH/dk_a is
    // the pos_matrix velocity operator: the diagonal kinetic 2 tpiba^2
    // (k+G)_a plus the separable projector derivative
    // (build_vkb/build_vkb_dk). The solved vector carries the complete
    // conduction-space position response and replaces the
    // empty-eigenvector-truncated r-matrix contraction.
    if (!wired() || hamilt_ == nullptr)
    {
        ModuleBase::timer::end("DFPT_PW", "solve_pos_resp");
        return;
    }
    const ModuleBase::Vector3<double> q_cart = data_.get_qvec(q_idx) * ucell_->G;
    const int nk = gs_psi_.get_nk();
    const int nbands = gs_psi_.get_nbands();
    const int lin_max = data_.get_max_iter();
    const double lin_thr = data_.get_conv_thr();

    for (int a = 0; a < 3; ++a)
    {
        std::vector<std::vector<std::vector<std::complex<double>>>> yvec(
            nk,
            std::vector<std::vector<std::complex<double>>>(nbands));
        for (int ik = 0; ik < nk; ++ik)
        {
            if (occ_kq_[ik].empty())
            {
                continue; // matches the displacement solve guard
            }
            if (last_q_ != q_idx || last_ik_ != ik)
            {
                hamilt_->set_context(q_cart, ik);
                last_q_ = q_idx;
                last_ik_ = ik;
            }
            pos_per_band_solve(ik, a, nbands, lin_max, lin_thr, yvec);
        }
        data_.set_pos_resp(a, yvec);
    }
    ModuleBase::timer::end("DFPT_PW", "solve_pos_resp");
}

void DFPT_PW::Impl::solve_efield_resp(int q_idx)
{
    ModuleBase::TITLE("DFPT_PW", "solve_efield_resp");
    ModuleBase::timer::start("DFPT_PW", "solve_efield_resp");
    // E-field SCF response (QE solve_e + dfpt_kernel form): the bare legs
    // Y^a stashed by solve_pos_resp are the field rhs base and the fixed
    // point adds the screened response potential of the mixed drho^E
    // exactly like solve_displacement. The converged dpsi^E,a feeds the
    // SCF dielectric tensor (DFPT_Q0::compute_eps) and the zstar_eu
    // cross-check probe (DFPT_ALEG).
    if (!wired() || hamilt_ == nullptr)
    {
        ModuleBase::timer::end("DFPT_PW", "solve_efield_resp");
        return;
    }
    const ModuleBase::Vector3<double> q_cart = data_.get_qvec(q_idx) * ucell_->G;
    const int npw = pw_rho_->npw;
    const int nk = gs_psi_.get_nk();
    const int nbands = gs_psi_.get_nbands();
    const int lin_max = data_.get_max_iter();
    const double lin_thr = data_.get_conv_thr();

    for (int a = 0; a < 3; ++a)
    {
        const std::vector<std::vector<std::vector<std::complex<double>>>> yr = data_.get_pos_resp(a);
        if (static_cast<int>(yr.size()) != nk)
        {
            continue; // bare legs not solved: no E response either
        }
        rho_.reset_mixing(q_idx);
        data_.set_drho_g(q_idx, 0, std::vector<std::complex<double>>(npw, std::complex<double>(0.0, 0.0)));
        bool converged = false;
        for (int iter = 0; iter < max_iter_ && !converged; ++iter)
        {
            // screened response potential of the mixed input density
            // (shared Impl::assemble_v_sc from dfpt_pw_solve.cpp)
            std::vector<std::complex<double>> v_sc_r;
            const std::vector<std::complex<double>> drho_in_g = data_.get_drho_g(q_idx, 0);
            assemble_v_sc(q_cart, drho_in_g, v_sc_r);
            for (int ik = 0; ik < nk; ++ik)
            {
                if (static_cast<int>(occ_kq_.size()) <= ik || occ_kq_[ik].empty())
                {
                    continue;
                }
                std::vector<std::vector<std::complex<double>>> dv_sc;
                pert_.apply_vr(q_idx, ik, v_sc_r, gs_psi_, q_cart, dv_sc);
                if (last_q_ != q_idx || last_ik_ != ik)
                {
                    hamilt_->set_context(q_cart, ik);
                    last_q_ = q_idx;
                    last_ik_ = ik;
                }
                efield_per_band_solve(ik, a, nbands, lin_max, lin_thr, yr, dv_sc);
            }
            rho_.compute_drho(gs_psi_, wg_, q_idx, data_);
            rho_.mix_drho(q_idx, data_);
            const double residual = rho_.get_residual(q_idx, data_);
            converged = (residual < conv_thr_);
            if (converged)
            {
                std::cout << "DFPT efield dir=" << a << " converged, residual=" << residual << " (iter=" << iter << ")"
                          << std::endl;
            }
        }
        stash_dpsi_efield(q_idx, a, nk, nbands);
    }
    ModuleBase::timer::end("DFPT_PW", "solve_efield_resp");
}

} // namespace ModuleDFPT
