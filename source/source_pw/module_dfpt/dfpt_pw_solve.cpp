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
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_base/tool_title.h"

#include <cmath>
#include <complex>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <sstream>
#include <iostream>
#include <vector>

namespace ModuleDFPT
{

namespace
{
/// relative defect ||out - in|| / ||out|| with convergence-safe guards
/// (shared definition with dfpt_rho.cpp; see the note there): a zero
/// output norm only reports convergence when the defect itself is zero,
/// and non-finite norms never report convergence
double safe_relative_defect(double dn2, double o2)
{
    if (!std::isfinite(dn2) || !std::isfinite(o2))
    {
        return std::numeric_limits<double>::max();
    }
    if (o2 <= 0.0)
    {
        return (dn2 <= 0.0) ? 0.0 : std::numeric_limits<double>::max();
    }
    return std::sqrt(dn2 / o2);
}
} // namespace

void DFPT_PW::Impl::assemble_v_sc(int q_idx,
                                  const ModuleBase::Vector3<double>& q_cart,
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
    // per-q context: analytic GGA kernels act with i(G + q) on the
    // physical carrier of the q-shifted amplitude
    xc_->set_context(q_cart);
    std::vector<std::complex<double>> a_r(nrxx);
    pw_rho_->recip2real(drho_in_g.data(), a_r.data());
    // nlcc: the first-order core density enters the XC channel only (the
    // core Coulomb is already in the ionic pseudopotential, so the Hartree
    // above correctly uses the valence density alone). The q index is
    // passed explicitly: the drhoc slot is per q and last_q_ can be stale
    // on the first assembly of a new q.
    const std::vector<std::complex<double>> drhoc_g = data_.get_drhoc_g(q_idx);
    if (static_cast<int>(drhoc_g.size()) == pw_rho_->npw)
    {
        std::vector<std::complex<double>> drhoc_r(nrxx);
        pw_rho_->recip2real(drhoc_g.data(), drhoc_r.data());
        for (int ir = 0; ir < nrxx; ++ir)
        {
            a_r[ir] += drhoc_r[ir];
        }
    }
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
    // stash the converged density response too: the two-pass accumulation
    // rebuilds f_xc.(drho_val + drhoc) per displacement after the working
    // drho slots have moved on; only needed when some species carries
    // a nonlinear core correction (skip the 3*nat*npw complex copies
    // for coreless runs)
    bool nlcc_any_solve = false;
    for (int it = 0; it < ucell_->ntype; ++it)
    {
        if (ucell_->atoms[it].ncpp.nlcc)
        {
            nlcc_any_solve = true;
            break;
        }
    }
    if (nlcc_any_solve)
    {
        const std::vector<std::complex<double>> drho_conv = data_.get_drho_g(q_idx, 0);
        data_.set_drho_disp(iat, idir, drho_conv);
    }
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
    // nlcc: q-shifted core-density derivative for the displaced atom
    // (fixed per displacement, consumed by assemble_v_sc's XC channel)
    rho_.compute_drhoc(*ucell_, iat, idir, q_idx, data_);
    rho_.reset_mixing(q_idx);
    data_.set_drho_g(q_idx, 0, std::vector<std::complex<double>>(npw, std::complex<double>(0.0, 0.0)));

    const int lin_max = data_.get_max_iter();
    const double lin_thr = data_.get_conv_thr();

    bool converged = false;
    double residual = 0.0;
    int n_iter_done = 0;
    std::vector<std::complex<double>> v_sc_r_last;
    // v_sc mode: the SCF mixes the screening potential instead of the
    // density. For a fixed linear kernel the density and potential
    // fixed-point maps have similar Jacobians (chi*K vs K*chi on the
    // responsive subspace), so scalar plain mixing keeps the same
    // spectrum: the variable change alone is not an inverse-kernel
    // preconditioner. What changes is the metric the Broyden least
    // squares works in (the K-weighted inner product), which is why
    // v_sc_broyden is the stabilized default while v_sc_plain tracks
    // the plain density mixing
    const bool vsc_mode = (mix_type_ == "v_sc_plain" || mix_type_ == "v_sc_broyden");
    std::vector<std::complex<double>> v_sc_in_r(nrxx, std::complex<double>(0.0, 0.0));
    // v_sc Broyden: a separate instance from the density Broyden because
    // the array length differs (nrxx real-space potential vs npw G-space
    // density); init and bind the Mixing_Data once per displacement cycle
    const bool vsc_broyden = (mix_type_ == "v_sc_broyden");
    if (vsc_broyden && !vsc_broyden_)
    {
        vsc_broyden_.reset(new Base_Mixing::Broyden_Mixing(broyden_ndim_));
    }
    if (vsc_broyden)
    {
        vsc_broyden_->reset();
        vsc_broyden_mdata_.reset(new Base_Mixing::Mixing_Data());
        vsc_broyden_->init_mixing_data(*vsc_broyden_mdata_, nrxx, sizeof(std::complex<double>));
    }
    for (int iter = 0; iter < max_iter_ && !converged; ++iter)
    {
        std::vector<std::complex<double>> v_sc_r;
        if (vsc_mode)
        {
            // use the mixed screening potential directly (the density is
            // only rebuilt from the converged potential for the 2n+1)
            v_sc_r = v_sc_in_r;
        }
        else
        {
            const std::vector<std::complex<double>> drho_in_g = data_.get_drho_g(q_idx, 0);
            assemble_v_sc(q_idx, q_cart, drho_in_g, v_sc_r);
        }
        if (dbg)
        {
            debug_iter_snapshot(iter, data_.get_drho_g(q_idx, 0), v_sc_r, npw, nrxx);
        }
        v_sc_r_last = vsc_mode ? v_sc_in_r : v_sc_r;

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
        if (vsc_mode)
        {
            // v_sc mode: rebuild v_sc from the OUTPUT density, then mix
            // the potential (not the density); the residual measures the
            // potential self-consistency
            std::vector<std::complex<double>> v_sc_out_r;
            const std::vector<std::complex<double>> drho_out_g = data_.get_drho_g(q_idx, 0);
            assemble_v_sc(q_idx, q_cart, drho_out_g, v_sc_out_r);
            if (vsc_broyden)
            {
                // Broyden on the real-space potential: push, cal_coef, mix
                auto screen = std::function<void(std::complex<double>*)>();
                const double beta = mix_beta_;
                auto mix_cb =
                    std::function<void(std::complex<double>*, const std::complex<double>*, const std::complex<double>*)>(
                        [nrxx, beta](std::complex<double>* out,
                                     const std::complex<double>* in,
                                     const std::complex<double>* F) {
                            for (int i = 0; i < nrxx; ++i)
                            {
                                out[i] = in[i] + beta * F[i];
                            }
                        });
                vsc_broyden_->push_data(*vsc_broyden_mdata_,
                                        v_sc_in_r.data(),
                                        v_sc_out_r.data(),
                                        screen,
                                        mix_cb,
                                        true);
                vsc_broyden_->cal_coef(
                    *vsc_broyden_mdata_,
                    std::function<double(std::complex<double>*, std::complex<double>*)>(
                        [nrxx](std::complex<double>* a, std::complex<double>* b) {
                            double s = 0.0;
                            for (int i = 0; i < nrxx; ++i)
                            {
                                s += (a[i] * std::conj(b[i])).real();
                            }
#ifdef __MPI
                            // pool-global inner product: each rank sums
                            // its FFT slab share only, and the Gram system
                            // must be identical on every rank
                            Parallel_Reduce::reduce_pool(s);
#endif
                            return s;
                        }));
                // residual from the unmixed in/out difference; the norms
                // are pool-global quantities (rank-local FFT slabs)
                double dn2 = 0.0;
                double o2 = 0.0;
                for (int ir = 0; ir < nrxx; ++ir)
                {
                    dn2 += std::norm(v_sc_out_r[ir] - v_sc_in_r[ir]);
                    o2 += std::norm(v_sc_out_r[ir]);
                }
#ifdef __MPI
                Parallel_Reduce::reduce_pool(dn2);
                Parallel_Reduce::reduce_pool(o2);
#endif
                residual = safe_relative_defect(dn2, o2);
                vsc_broyden_->mix_data(*vsc_broyden_mdata_, v_sc_in_r.data());
            }
            else
            {
                // plain beta mix on the real-space potential; the norms are
                // pool-global quantities (rank-local FFT slabs)
                double dn2 = 0.0;
                double o2 = 0.0;
                for (int ir = 0; ir < nrxx; ++ir)
                {
                    const std::complex<double> delta = v_sc_out_r[ir] - v_sc_in_r[ir];
                    dn2 += std::norm(delta);
                    o2 += std::norm(v_sc_out_r[ir]);
                    v_sc_in_r[ir] += mix_beta_ * delta;
                }
#ifdef __MPI
                Parallel_Reduce::reduce_pool(dn2);
                Parallel_Reduce::reduce_pool(o2);
#endif
                residual = safe_relative_defect(dn2, o2);
            }
            // the 2n+1 accumulation path reads the density from data_:
            // unlike the density-mixing branch, the v_sc loop never calls
            // mix_drho, so data_ still holds the last compute_drho output
            // (drho_out_g was read from it), which at convergence is the
            // density consistent with the final screening potential
        }
        else
        {
            rho_.mix_drho(q_idx, data_);
            residual = rho_.get_residual(q_idx, data_);
        }
        // default iteration trace (GS-SCF style); the DFPT_DEBUG snapshots
        // above stay opt-in diagnostics. Formatting goes through a local
        // buffer so the stream flags/precision are never disturbed for
        // the surrounding output (#8164 review: setprecision pollution
        // degraded the DBG snapshots of the next iteration)
        {
            std::ostringstream trace;
            trace << "  DFPT mix iter=" << std::setw(3) << std::right << iter + 1 << "/" << max_iter_ << " residual="
                  << std::setprecision(3) << std::scientific << residual << " (thr " << std::setprecision(1)
                  << std::scientific << conv_thr_ << ")" << std::endl;
            std::cout << trace.str();
        }
        if (dbg)
        {
            std::cout << "DBG iter=" << iter << " residual=" << residual << " conv_thr=" << conv_thr_ << std::endl;
        }
        converged = (residual < conv_thr_);
        ++n_iter_done;
    }
    // convergence summary in the efield-solve style; names the state
    // explicitly instead of failing silently at max_iter; formatted in a
    // local buffer for the same stream-state reason
    {
        std::ostringstream summary;
        summary << "DFPT disp q=(" << q_frac.x << "," << q_frac.y << "," << q_frac.z << ") atom=" << iat + 1 << "("
                << ucell_->atoms[ucell_->iat2it[iat]].label << ") dir=" << static_cast<char>('x' + idir) << " "
                << (converged ? "converged" : "NOT converged (max_iter reached)") << ", residual=" << std::scientific
                << std::setprecision(6) << residual << " (iter=" << n_iter_done << ")" << std::endl;
        std::cout << summary.str();
    }
    if (!converged)
    {
        std::ostringstream oss;
        oss << "the displacement response (atom " << iat + 1 << " dir " << static_cast<char>('x' + idir)
            << ") did not converge within dfpt_max_iter=" << max_iter_ << ": final residual " << std::scientific
            << std::setprecision(3) << residual << " vs threshold " << std::setprecision(1) << conv_thr_
            << std::defaultfloat << ". Phonon frequencies involving this displacement are unreliable";
        ModuleBase::WARNING("DFPT_PW::solve_displacement", oss.str());
    }
    stash_converged_disp_response(q_idx, iat, idir, v_sc_r_last, nk, nbands);

    ModuleBase::timer::end("DFPT_PW", "solve_displacement");
    return residual;
}

} // namespace ModuleDFPT
