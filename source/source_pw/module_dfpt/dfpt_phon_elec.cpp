// Electronic (2n+1) contribution of DFPT_Phon::accumulate_electron,
// split out of dfpt_phon.cpp: the converged-dpsi stash/restore, the
// bare-potential cross sum and the same-atom second-order term. All
// formulas are moved verbatim from the original body.

#include "dfpt_phon.h"

#include "dfpt_kq_basis.h"
#include "dfpt_pert.h"
#include "dfpt_pw_data.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>

#include "source_base/timer.h"
#include "source_base/tool_title.h"

namespace ModuleDFPT
{

namespace
{

/// dpsi layout of one displacement: [k][band][G]
using DpsiDisp = std::vector<std::vector<std::vector<std::complex<double>>>>;

/// stash the converged dpsi of one displacement: prefer the per-displacement
/// store of the two-pass flow; fall back to the working slots for the legacy
/// interleaved call order
DpsiDisp stash_dpsi(DFPT_PW_Data& data, int q_idx, int atom_idx, int dir, int nk, int nbands)
{
    DpsiDisp dpsib = data.get_dpsi_disp(atom_idx, dir);
    if (dpsib.empty() || static_cast<int>(dpsib.size()) < nk
        || (nk > 0 && static_cast<int>(dpsib[0].size()) < nbands))
    {
        dpsib.assign(nk, std::vector<std::vector<std::complex<double>>>(nbands));
        for (int ik = 0; ik < nk; ++ik)
        {
            for (int ib = 0; ib < nbands; ++ib)
            {
                dpsib[ik][ib] = data.get_dpsi(q_idx, ik, ib);
            }
        }
    }
    return dpsib;
}

/// restore the stashed converged dpsi after the accumulation touched the
/// working slots (apply_dv reuses the slot of every k/band)
void restore_dpsi(DFPT_PW_Data& data, int q_idx, int nk, int nbands, const DpsiDisp& dpsib)
{
    for (int ik = 0; ik < nk; ++ik)
    {
        for (int ib = 0; ib < nbands; ++ib)
        {
            if (!dpsib[ik][ib].empty())
            {
                data.set_dpsi(q_idx, ik, ib, dpsib[ik][ib]);
            }
        }
    }
}

/// bare-potential cross sum X_ba = sum_kn wg <dpsi^b | dV^a_ext | psi>;
/// dV^a_ext must have been built (build_dv) before this call
std::complex<double> cross_sum(DFPT_Pert& pert,
                               DFPT_PW_Data& data,
                               int q_idx,
                               const DpsiDisp& dpsib,
                               const psi::Psi<std::complex<double>>& psi,
                               const ModuleBase::matrix& wg)
{
    const int nk = psi.get_nk();
    const int nbands = psi.get_nbands();
    std::complex<double> cross(0.0, 0.0);
    for (int ik = 0; ik < nk; ++ik)
    {
        pert.apply_dv(q_idx, ik, psi, data);
        for (int ib = 0; ib < nbands; ++ib)
        {
            if (!dfpt_band_occupied(wg, ik, ib))
            {
                continue;
            }
            const std::vector<std::complex<double>> rhs = data.get_dpsi(q_idx, ik, ib);
            const std::vector<std::complex<double>>& sol = dpsib[ik][ib];
            if (rhs.size() != sol.size() || sol.empty())
            {
                continue;
            }
            std::complex<double> dot(0.0, 0.0);
            for (size_t i = 0; i < sol.size(); ++i)
            {
                dot += std::conj(sol[i]) * rhs[i];
            }
            cross += wg(ik, ib) * dot;
        }
    }
    return cross;
}

/// scatter the band-ib q-carrier chi from the k+q_eff G-shell onto the rho
/// grid and transform to real space; a carrier of the wrong shape (chi not
/// stored for this k/band) leaves x_r zero, matching apply_d2vnl
void scatter_chi_to_r(const ModulePW::PW_Basis& pw_rho,
                      const std::vector<std::vector<std::complex<double>>>& chi,
                      int nbands,
                      int ib,
                      int npwk_kq,
                      const DFPT_KQ_Basis& kq,
                      std::vector<std::complex<double>>& x_r)
{
    std::vector<std::complex<double>> x_recip(pw_rho.npw, std::complex<double>(0.0, 0.0));
    if (static_cast<int>(chi.size()) == nbands && static_cast<int>(chi[ib].size()) == npwk_kq)
    {
        for (int igl = 0; igl < npwk_kq; ++igl)
        {
            const int ig_rho = kq.get_ig_rho(igl);
            if (ig_rho >= 0)
            {
                x_recip[ig_rho] = chi[ib][igl];
            }
        }
        pw_rho.recip2real(x_recip.data(), x_r.data());
    }
    else
    {
        std::fill(x_r.begin(), x_r.end(), std::complex<double>(0.0, 0.0));
    }
}

} // namespace

void DFPT_Phon::accumulate_electron(int q_idx,
                                    int atom_idx,
                                    int dir,
                                    const psi::Psi<std::complex<double>>& psi,
                                    const ModuleBase::matrix& wg,
                                    DFPT_PW_Data& data)
{
    ModuleBase::TITLE("DFPT_Phon", "accumulate_electron");
    ModuleBase::timer::start("DFPT_Phon", "accumulate_electron");
    if (pert_ == nullptr || pw_rho_ == nullptr || ucell_ == nullptr)
    {
        ModuleBase::timer::end("DFPT_Phon", "accumulate_electron");
        return;
    }
    const int nat = ucell_->nat;
    const int nat3 = 3 * nat;
    if (accum_q_ != q_idx || dynmat_accum_.nr != nat3)
    {
        dynmat_accum_ = ModuleBase::ComplexMatrix(nat3, nat3, true);
        accum_q_ = q_idx;
    }
    const int rowb = 3 * atom_idx + dir;

    // stash the converged dpsi of this displacement (apply_dv reuses the slot)
    const DpsiDisp dpsib = stash_dpsi(data, q_idx, atom_idx, dir, psi.get_nk(), psi.get_nbands());

    for (int iat = 0; iat < nat; ++iat)
    {
        for (int idir = 0; idir < 3; ++idir)
        {
            const int cola = 3 * iat + idir;
            // ---- term 2 <dpsi^b | dV^a_ext | psi> over all k,n ----
            // Hermitian (2n+1) accumulation: the row element gets X_ba and
            // the transposed element gets conj(X_ba); the self-consistent
            // response of dpsi^b already contains the screening, and the
            // Hartree-xc kernel quadratic term cancels the <dpsi|dV_sc|psi>
            // cross terms by the variational identity, so only the bare
            // external perturbation appears here
            pert_->build_dv(q_idx, iat, idir, data);
            const std::complex<double> cross = cross_sum(*pert_, data, q_idx, dpsib, psi, wg);
            const double mass_norm
                = std::sqrt(ucell_->atoms[ucell_->iat2it[atom_idx]].mass * ucell_->atoms[ucell_->iat2it[iat]].mass);
            dynmat_accum_(rowb, cola) += cross / mass_norm;
            dynmat_accum_(cola, rowb) += std::conj(cross) / mass_norm;

            // ---- same-atom anharmonic term <psi | d2_ab V_ext | psi> ----
            if (iat == atom_idx)
            {
                accum_d2_same_atom(atom_idx, dir, iat, idir, rowb, psi, wg);
            }
        }
    }

    // restore the converged dpsi of this displacement
    restore_dpsi(data, q_idx, psi.get_nk(), psi.get_nbands(), dpsib);
    ModuleBase::timer::end("DFPT_Phon", "accumulate_electron");
}

void DFPT_Phon::accum_d2_same_atom(int atom_idx,
                                   int dir,
                                   int iat,
                                   int idir,
                                   int rowb,
                                   const psi::Psi<std::complex<double>>& psi,
                                   const ModuleBase::matrix& wg)
{
    const int cola = 3 * iat + idir;
    if (cola < rowb)
    {
        // only the upper triangle of the same-atom block is accumulated
        // here; the Hermitian partner is added below from conj(d2sum)
        return;
    }
    const int nk = psi.get_nk();
    const int nbands = psi.get_nbands();
    // QE ground truth (dynmat_us.f90 + phq_init.f90): the mixed
    // (+q, -q) second-order potential of the local part is
    // -Omega tpiba^2 G_a G_b vloc(|G|) Re[rho(G) e^{-iG tau_s}]
    // (integer G, no q), and the KB nonlocal part is the same-atom
    // block deff[gammap*becp1 + becp1*gammap + alphap_a*alphap_b +
    // alphap_b*alphap_a] with becp1/alphap/gammap all built from
    // vkb_k and (k+G) factors (integer G, no q). The (+q,-q)
    // dressings collapse to an integer-G carrier for every q, so
    // this term is q-independent and must never be gated on 2q
    // commensurability (the old gate silently dropped it for
    // 2q not reciprocal, e.g. q=(0.25,0,0), and produced
    // imaginary phonon branches).
    const ModuleBase::Vector3<double> q_eff_cart(0.0, 0.0, 0.0);
    std::vector<std::complex<double>> dv2_r;
    pert_->d2vloc_r(atom_idx, idir, dir, dv2_r);
    if (static_cast<int>(dv2_r.size()) != pw_rho_->nrxx)
    {
        dv2_r.assign(pw_rho_->nrxx, std::complex<double>(0.0, 0.0));
    }
    std::vector<std::vector<std::complex<double>>> chi;
    std::complex<double> d2sum(0.0, 0.0);
    std::vector<std::complex<double>> u_r(pw_rho_->nrxx);
    std::vector<std::complex<double>> x_r(pw_rho_->nrxx);
    for (int ik = 0; ik < nk; ++ik)
    {
        pert_->apply_d2vnl(atom_idx, idir, dir, q_eff_cart, psi, ik, chi);
        // k+q_eff scatter map for this k (must match apply_d2vnl)
        DFPT_KQ_Basis kq;
        kq.init(pert_->get_pw_wfc(), pert_->get_pw_rho(), q_eff_cart, ik);
        const int npwk_kq = kq.get_npwk();
        for (int ib = 0; ib < nbands; ++ib)
        {
            if (!dfpt_band_occupied(wg, ik, ib))
            {
                continue;
            }
            pert_->get_pw_wfc()->recip2real(&psi(ik, ib, 0), u_r.data(), ik);
            scatter_chi_to_r(*pw_rho_, chi, nbands, ib, npwk_kq, kq, x_r);
            std::complex<double> expect(0.0, 0.0);
            for (int ir = 0; ir < pw_rho_->nrxx; ++ir)
            {
                expect += std::conj(u_r[ir]) * u_r[ir] * dv2_r[ir] + std::conj(u_r[ir]) * x_r[ir];
            }
            d2sum += wg(ik, ib) * expect / static_cast<double>(pw_rho_->nxyz);
        }
    }
    const double inv_m = 1.0 / ucell_->atoms[ucell_->iat2it[atom_idx]].mass;
    dynmat_accum_(rowb, cola) += d2sum * inv_m;
    if (cola != rowb)
    {
        dynmat_accum_(cola, rowb) += std::conj(d2sum) * inv_m;
    }
}

} // namespace ModuleDFPT
