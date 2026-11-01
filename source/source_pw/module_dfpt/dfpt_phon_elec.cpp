// Electronic (2n+1) contribution of DFPT_Phon::accumulate_electron,
// split out of dfpt_phon.cpp: the converged-dpsi stash/restore, the
// bare-potential cross sum and the same-atom second-order term. All
// formulas are moved verbatim from the original body.

#include "dfpt_phon.h"

#include "dfpt_kq_basis.h"
#include "dfpt_pert.h"
#include "dfpt_pw_data.h"
#include "dfpt_rho.h"
#include "source_base/math_integral.h"
#include "source_base/math_sphbes.h"
#include "source_base/parallel_reduce.h"

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
                                     DFPT_PW_Data& data,
                                     const std::vector<std::complex<double>>& dvxc_r)
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

    // ---- nlcc direct terms (QE dynmat_nlcc + dynmatcc) ----
    accum_nlcc_kernel(atom_idx, dir, q_idx, data, dvxc_r);
    if (dir == 0)
    {
        // the diagonal block is direction-pair complete: add it once per
        // atom (every accumulate_electron call would triple it)
        accum_nlcc_diag(atom_idx, q_idx);
    }

    // restore the converged dpsi of this displacement
    restore_dpsi(data, q_idx, psi.get_nk(), psi.get_nbands(), dpsib);
    ModuleBase::timer::end("DFPT_Phon", "accumulate_electron");
}

void DFPT_Phon::accum_nlcc_kernel(int atom_idx,
                                   int dir,
                                   int q_idx,
                                   DFPT_PW_Data& data,
                                   const std::vector<std::complex<double>>& dvxc_r)
{
    ModuleBase::TITLE("DFPT_Phon", "accum_nlcc_kernel");
    ModuleBase::timer::start("DFPT_Phon", "accum_nlcc_kernel");
    const int nrxx = pw_rho_->nrxx;
    if (rho_ == nullptr || ucell_ == nullptr || static_cast<int>(dvxc_r.size()) != nrxx)
    {
        ModuleBase::timer::end("DFPT_Phon", "accum_nlcc_kernel");
        return;
    }
    // any nlcc atom at all?
    bool nlcc_any = false;
    for (int it = 0; it < ucell_->ntype; ++it)
    {
        if (ucell_->atoms[it].ncpp.nlcc)
        {
            nlcc_any = true;
            break;
        }
    }
    if (!nlcc_any)
    {
        ModuleBase::timer::end("DFPT_Phon", "accum_nlcc_kernel");
        return;
    }
    const int nat = ucell_->nat;
    const ModuleBase::Vector3<double> q_cart = data.get_qvec(q_idx) * ucell_->G;
    // pool-global normalization: nrxx is the rank-local FFT slab size,
    // so the raw sum needs a pool reduction (the mixing PR's convention)
    const double fac = ucell_->omega / static_cast<double>(pw_rho_->nxyz);
    const int rowb = 3 * atom_idx + dir;
    const double mass_row = ucell_->atoms[ucell_->iat2it[atom_idx]].mass;
    // QE dynmat_nlcc semantics, pinned by the element-level comparison
    // against QE dyn.G (bn_cubic): the raw matrix carries
    //   X_ab = (Omega/nrxx) sum_r dvxc^a(r) conj(drhoc^b(r)),
    // dvxc^a = f_xc.(drho_val^a + drhoc^a), once per element, and the
    // effective Hermitian value is the symmetrized
    //   S_ab = (1/2)(X_ab + conj(X_ba)).
    // X_ab is two-sided only when BOTH displacements carry core response;
    // for the one-sided cross blocks (e.g. B vs N: X_BN = 0 because N has
    // no nlcc) this reduces to (1/2) X_NB. The column loop must therefore
    // cover non-nlcc atoms too: their X_ab vanishes but the transposed
    // piece conj(X_ba) = sum_r dvxc^b conj(drhoc^a) -- the column's
    // valence-only dvaux against the ROW's core -- is exactly what the
    // (nlcc-row, non-nlcc-column) element needs.
    std::vector<std::complex<double>> drhoc_g;
    std::vector<std::complex<double>> drhoc_r(nrxx);
    // row's own core response in real space (for the conj(X_ba) piece;
    // empty when the row atom carries no nlcc)
    std::vector<std::complex<double>> drhoc_a_r;
    rho_->compute_drhoc_g(*ucell_, atom_idx, dir, q_cart, drhoc_g);
    const bool row_has_core = static_cast<int>(drhoc_g.size()) == pw_rho_->npw;
    if (row_has_core)
    {
        drhoc_a_r.resize(nrxx);
        pw_rho_->recip2real(drhoc_g.data(), drhoc_a_r.data());
    }
    const bool have_table = dvxc_core_ != nullptr;
    for (int iat2 = 0; iat2 < nat; ++iat2)
    {
        const bool col_has_core = ucell_->atoms[ucell_->iat2it[iat2]].ncpp.nlcc;
        for (int idir2 = 0; idir2 < 3; ++idir2)
        {
            // X_ab: the row's dvaux against the column's core response.
            // Conjugation convention: the ROW response carries the
            // conjugate (sum(conj(row)*col), matching the wavefunction
            // and Ewald terms; at Gamma with real amplitudes this is a
            // no-op, at nonzero q it fixes the imaginary cross terms)
            std::complex<double> x_ab(0.0, 0.0);
            if (col_has_core)
            {
                rho_->compute_drhoc_g(*ucell_, iat2, idir2, q_cart, drhoc_g);
                if (static_cast<int>(drhoc_g.size()) != pw_rho_->npw)
                {
                    continue;
                }
                pw_rho_->recip2real(drhoc_g.data(), drhoc_r.data());
                for (int ir = 0; ir < nrxx; ++ir)
                {
                    x_ab += std::conj(dvxc_r[ir]) * drhoc_r[ir];
                }
#ifdef __MPI
                Parallel_Reduce::reduce_pool(x_ab);
#endif
                x_ab *= fac;
            }
            // X_ba: the column's dvaux (per-q table) against the row's
            // core response; same conjugation swap
            std::complex<double> x_ba(0.0, 0.0);
            const int cola = 3 * iat2 + idir2;
            if (row_has_core && have_table && cola < static_cast<int>(dvxc_core_->size())
                && static_cast<int>((*dvxc_core_)[cola].size()) == nrxx)
            {
                const std::vector<std::complex<double>>& dvxc_b = (*dvxc_core_)[cola];
                for (int ir = 0; ir < nrxx; ++ir)
                {
                    x_ba += std::conj(dvxc_b[ir]) * drhoc_a_r[ir];
                }
#ifdef __MPI
                Parallel_Reduce::reduce_pool(x_ba);
#endif
                x_ba *= fac;
            }
            if (x_ab == std::complex<double>(0.0, 0.0) && x_ba == std::complex<double>(0.0, 0.0))
            {
                continue;
            }
            const std::complex<double> s_ab = 0.5 * (x_ab + std::conj(x_ba));
            const double mass_norm = std::sqrt(mass_row * ucell_->atoms[ucell_->iat2it[iat2]].mass);
            // single-ordering fill: the transposed element receives its
            // own S_ba = conj(S_ab) from the column displacement's pass
            dynmat_accum_(rowb, cola) += s_ab / mass_norm;
        }
    }
    ModuleBase::timer::end("DFPT_Phon", "accum_nlcc_kernel");
}

void DFPT_Phon::accum_nlcc_diag(int atom_idx, int q_idx)
{
    ModuleBase::TITLE("DFPT_Phon", "accum_nlcc_diag");
    ModuleBase::timer::start("DFPT_Phon", "accum_nlcc_diag");
    const int npw = pw_rho_ ? pw_rho_->npw : 0;
    const int nrxx = pw_rho_ ? pw_rho_->nrxx : 0;
    if (ucell_ == nullptr || pw_rho_ == nullptr || vxc_r_ == nullptr
        || static_cast<int>(vxc_r_->size()) != nrxx || npw == 0 || accum_q_ != q_idx)
    {
        ModuleBase::timer::end("DFPT_Phon", "accum_nlcc_diag");
        return;
    }
    const int it = ucell_->iat2it[atom_idx];
    if (!ucell_->atoms[it].ncpp.nlcc)
    {
        ModuleBase::timer::end("DFPT_Phon", "accum_nlcc_diag");
        return;
    }
    // QE dynmatcc ground truth, evaluated as a real-space integral:
    //   D_ij += Int dr d2rho_c/(d tau_i d tau_j) (r) * v_xc^0(r)
    // with the second-order shifted-core amplitude
    //   d2rho_c(G) = -tpiba^2 G_i G_j rhoc~(|G|) exp(-i G tau_a)
    // (a real field: W(-G) = conj(W(G))), so the contraction runs through
    // the same recip2real machinery as every other real-space field and
    // never touches the half-space G-sum convention directly.
    const pseudo& p = ucell_->atoms[it].ncpp;
    const int mesh = p.msh;
    const ModuleBase::Vector3<double> tau_cart
        = ucell_->atoms[it].tau[ucell_->iat2ia[atom_idx]] * ucell_->lat0;
    std::vector<double> aux(mesh);
    std::vector<std::complex<double>> w_g(npw);
    std::vector<std::complex<double>> w_r(nrxx);
    double d3[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    const double fac = ucell_->omega / static_cast<double>(pw_rho_->nxyz);
    for (int i = 0; i < 3; ++i)
    {
        for (int j = i; j < 3; ++j)
        {
            std::fill(w_g.begin(), w_g.end(), std::complex<double>(0.0, 0.0));
            for (int ig = 0; ig < npw; ++ig)
            {
                const ModuleBase::Vector3<double>& gv = pw_rho_->gcar[ig];
                const double g2 = gv * gv;
                // the G = 0 shell contributes nothing (G_i G_j factor)
                if (g2 < 1.0e-12)
                {
                    continue;
                }
                const double gphys = std::sqrt(g2) * ucell_->tpiba;
                ModuleBase::Sphbes::Spherical_Bessel(mesh, p.r.data(), gphys, 0, aux.data());
                double rhocg = 0.0;
                for (int ir = 0; ir < mesh; ++ir)
                {
                    aux[ir] = p.r[ir] * p.r[ir] * p.rho_atc[ir] * aux[ir];
                }
                ModuleBase::Integral::Simpson_Integral(mesh, aux.data(), p.rab.data(), rhocg);
                rhocg *= ModuleBase::FOUR_PI / ucell_->omega;

                const double phase = ModuleBase::TWO_PI * (gv * tau_cart) / ucell_->lat0;
                // AGENTS rule 14: named locals instead of computed ctor args
                const double gtau_re = std::cos(phase);
                const double gtau_im = -std::sin(phase);
                const std::complex<double> gtau(gtau_re, gtau_im);
                w_g[ig] = -ucell_->tpiba * ucell_->tpiba * gv[i] * gv[j] * rhocg * gtau;
            }
            pw_rho_->recip2real(w_g.data(), w_r.data());
            double s = 0.0;
            for (int ir = 0; ir < nrxx; ++ir)
            {
                s += w_r[ir].real() * (*vxc_r_)[ir];
            }
#ifdef __MPI
            Parallel_Reduce::reduce_pool(s);
#endif
            d3[i][j] = fac * s;
            if (j != i)
            {
                d3[j][i] = d3[i][j];
            }
        }
    }
    const double mass = ucell_->atoms[it].mass;
    const int base = 3 * atom_idx;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            dynmat_accum_(base + i, base + j) += d3[i][j] / mass;
        }
    }
    ModuleBase::timer::end("DFPT_Phon", "accum_nlcc_diag");
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
