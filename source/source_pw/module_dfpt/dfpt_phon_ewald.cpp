// Ewald ion-ion part of DFPT_Phon::ion_ion, rewritten following the
// structure of QE 6.8 PHonon/PH/d2ionq.f90:
//   1. reciprocal-space sum over |G+q| kernels for ALL atom pairs
//      (including na == nb, where the phase is trivially 1), so the
//      diagonal blocks need no separate self machinery;
//   2. a Gamma-kernel subtraction on the diagonal blocks only, with the
//      structure-factor phases summed over ALL neighbor atoms (the
//      acoustic-sum-rule-preserving bookkeeping);
//   3. a real-space erfc-screened pair sum over the neighbor shells of the
//      shared H_Ewald_pw::rgen generator (same routine as the ground-state
//      energy/force/stress Ewald sums), with the full complex phase
//      e^{i 2pi q.L} on the pair blocks and the phase-free negative
//      counterpart accumulated on the row atom's diagonal in the same loop;
//      the integer triplet of each shell translation is recovered exactly
//      as n = R . G^T (a_i . g_j = delta_ij convention).
// The in-memory matrix keeps the ABACUS convention D = Phi/sqrt(Ma Mb)
// (Ry/bohr^2/amu): pair coefficients carry 1/sqrt(Ma Mb), diagonal-block
// coefficients carry 1/Ma.

#include "dfpt_phon.h"

#include "source_base/constants.h"
#include "source_base/parallel_reduce.h"
#include "source_base/tool_quit.h"
#include "source_base/truncated_func.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/unitcell.h"
#include "source_hamilt/module_ewald/h_ewald_pw.h"

#include <cmath>
#include <complex>
#include <iomanip>
#include <vector>

#include "source_base/timer.h"
#include "source_base/tool_title.h"

namespace ModuleDFPT
{

namespace
{

/// signed frequencies helper is not needed here; the radial kernels below
/// mirror QE d2ionq.f90: the second derivative of erfc(sqrt(a) r)/r splits
/// into a traceless part d2f (contracted with r_i r_j) and an isotropic
/// part df (contracted with delta_ij)
void ewald_radial(
    double ar,
    double rr,
    double &d2f,
    double &df)
{
    const double er = ModuleBase::truncated_erfc(ar);
    const double ex = std::exp(-ar * ar);
    const double sq8tpi = std::sqrt(8.0 / ModuleBase::TWO_PI);
    d2f = (3.0 * er + sq8tpi * ar * (3.0 + 2.0 * ar * ar) * ex) / (rr * rr * rr * rr * rr);
    df = (-er - sq8tpi * ar * ex) / (rr * rr * rr);
}

} // namespace

/// screening alpha of the Ewald split, searched so that the G-sum tail is
/// bounded inside the rho grid. Multiplicative shrinkage (the historical
/// schedule) can continue below 0.1, where the fixed-decrement variant
/// crosses zero and aborts even though a solution exists: for total
/// charge 8 and ggecut*tpiba2 = 8, alpha = 0.08 satisfies the 1e-9 tail
/// while alpha = 0.10 misses it.
double ewald_alpha_search(double charge, double ggecut, double tpiba2)
{
    const double alpha_init = 1.1;        ///< empirical parameter: initial Ewald screening-alpha guess
    const double alpha_shrink = 0.9;      ///< empirical parameter: per-iteration alpha shrink factor
    const double alpha_min = 1.0e-4;      ///< empirical parameter: alpha floor of the search
    const double ewald_tail_tol = 1.0e-9; ///< empirical parameter: accepted G-sum tail bound
    // ggecut counts |G_max|^2 in 1/lat0^2 units (pw_basis.h), so the bohr^2
    // cutoff is ggecut * tpiba2
    double alpha = alpha_init;
    double upperbound = 0.0;
    do
    {
        alpha *= alpha_shrink;
        if (alpha < alpha_min)
        {
            ModuleBase::WARNING_QUIT("DFPT_Phon::ion_ion", "Can't find optimal Ewald alpha.");
        }
        upperbound = 2.0 * charge * charge * std::sqrt(2.0 * alpha / ModuleBase::TWO_PI) *
                     ModuleBase::truncated_erfc(std::sqrt(ggecut * tpiba2 / 4.0 / alpha));
    } while (upperbound > ewald_tail_tol);
    return alpha;
}

void DFPT_Phon::ion_ion(const ModuleBase::Vector3<double> &q_frac, ModuleBase::ComplexMatrix &dyn)
{
    ModuleBase::TITLE("DFPT_Phon", "ion_ion");
    ModuleBase::timer::start("DFPT_Phon", "ion_ion");
    if (ucell_ == nullptr || pw_rho_ == nullptr)
    {
        ModuleBase::timer::end("DFPT_Phon", "ion_ion");
        return;
    }
    const int nat = ucell_->nat;
    const int nat3 = 3 * nat;

    // total ionic charge
    double charge = 0.0;
    for (int it = 0; it < ucell_->ntype; ++it)
    {
        charge += ucell_->atoms[it].na * ucell_->atoms[it].ncpp.zv;
    }

    // choose the screening alpha so the G-sum tail is converged inside the
    // rho grid; derive the real-space cutoff (QE: terms up to erfc(5) are
    // counted)
    ewald_alpha_ = ewald_alpha_search(charge, pw_rho_->ggecut, ucell_->tpiba2);
    if (ewald_alpha_override_ > 0.0)
    {
        ewald_alpha_ = ewald_alpha_override_;
        ewald_alpha_override_ = 0.0; // one-shot pin, production never touches it
    }
    const double ewald_rcut_factor = 5.0; ///< empirical parameter: real-space cutoff in 1/sqrt(alpha)
    ewald_rcut_ = ewald_rcut_factor / std::sqrt(ewald_alpha_);
    const double rmax_lat0 = ewald_rcut_ / ucell_->lat0;

    const ModuleBase::Vector3<double> q_cart = q_frac * ucell_->G;

    ModuleBase::ComplexMatrix dy1(nat3, nat3, true);
    ModuleBase::ComplexMatrix dy2(nat3, nat3, true);
    ModuleBase::ComplexMatrix dy3(nat3, nat3, true);

    // ---------------- reciprocal-space sums ----------------
    // dy1: |G+q| kernel, all atom pairs including na == nb (the phase
    // e^{i 2pi (G+q).(tau_a - tau_b)} is 1 on the diagonal)
    // dy2: Gamma kernel (G only, G != 0), cos phases summed over ALL
    // neighbor atoms, written only on the row atom's diagonal block
    const double w2_floor = 1.0e-12; ///< empirical parameter: |G+q|^2 zero-shell guard (1/lat0^2)
    const double g2_floor = 1.0e-12; ///< empirical parameter: |G|^2 zero-shell guard (1/lat0^2)
    for (int ig = 0; ig < pw_rho_->npw; ++ig)
    {
        const ModuleBase::Vector3<double> &gcart = pw_rho_->gcar[ig];
        const ModuleBase::Vector3<double> w = gcart + q_cart;
        const double w2 = w * w;
        const double g2 = gcart * gcart;
        double facq = 0.0;
        if (w2 > w2_floor)
        {
            const double w2_bohr = w2 * ucell_->tpiba2;
            facq =
                -ModuleBase::e2 * ModuleBase::FOUR_PI / ucell_->omega * std::exp(-w2_bohr / (4.0 * ewald_alpha_)) / w2;
        }
        double fac = 0.0;
        if (g2 > g2_floor)
        {
            const double g2_bohr = g2 * ucell_->tpiba2;
            fac =
                -ModuleBase::e2 * ModuleBase::FOUR_PI / ucell_->omega * std::exp(-g2_bohr / (4.0 * ewald_alpha_)) / g2;
        }
        for (int na = 0; na < nat; ++na)
        {
            const int ita = ucell_->iat2it[na];
            const int iia = ucell_->iat2ia[na];
            const double za = ucell_->atoms[ita].ncpp.zv;
            const double ma = ucell_->atoms[ita].mass;
            const ModuleBase::Vector3<double> &ta = ucell_->atoms[ita].tau[iia];
            for (int nb = 0; nb < nat; ++nb)
            {
                const int itb = ucell_->iat2it[nb];
                const int iib = ucell_->iat2ia[nb];
                const double zb = ucell_->atoms[itb].ncpp.zv;
                const double mb = ucell_->atoms[itb].mass;
                const ModuleBase::Vector3<double> &tb = ucell_->atoms[itb].tau[iib];
                // dy1: all pairs, mass-scaled pair coefficient 1/sqrt(Ma Mb)
                if (facq != 0.0)
                {
                    const double argq = ModuleBase::TWO_PI * (w * (ta - tb));
                    const std::complex<double> phase(std::cos(argq), std::sin(argq));
                    const std::complex<double> facg = facq * za * zb / std::sqrt(ma * mb) * phase;
                    for (int da = 0; da < 3; ++da)
                    {
                        for (int db = 0; db < 3; ++db)
                        {
                            dy1(3 * na + da, 3 * nb + db) += facg * w[da] * w[db];
                        }
                    }
                }
                // dy2: the Gamma-phase kernel on the row atom's diagonal
                // block, summed over ALL neighbor atoms (nb == na included
                // with cos 0 = 1); the diagonal-block mass scaling is 1/Ma
                if (fac != 0.0)
                {
                    const double arg = ModuleBase::TWO_PI * (gcart * (ta - tb));
                    const double fnat_coeff = fac * za * zb * std::cos(arg) / ma;
                    for (int da = 0; da < 3; ++da)
                    {
                        for (int db = 0; db < 3; ++db)
                        {
                            dy2(3 * na + da, 3 * na + db) += fnat_coeff * gcart[da] * gcart[db];
                        }
                    }
                }
            }
        }
    }

    // ---------------- real-space sum ----------------
    // shells R = n . latvec with |R - dtau| <= rmax carry the full complex
    // phase e^{i 2pi q.n} on the pair block; the phase-free negative
    // counterpart lands on the row atom's diagonal in the same loop
    // (QE d2ionq: only the rank owning G = 0 enters this part).
    // The neighbor shells come from the shared H_Ewald_pw::rgen generator
    // (same routine as the ground-state energy/force/stress Ewald sums):
    // rgen returns r = R - dtau in lat0 units, sorted by increasing length.
    const bool has_g0 = (pw_rho_->ig_gge0 >= 0);
    if (has_g0)
    {
        const int mxr = H_Ewald_pw::estimate_mxr(rmax_lat0, ucell_->G);
        std::vector<ModuleBase::Vector3<double>> r(mxr);
        std::vector<double> r2(mxr);
        std::vector<int> irr(mxr);
        int nrm = 0;
        const ModuleBase::Matrix3 &G = ucell_->G;
        for (int na = 0; na < nat; ++na)
        {
            const int ita = ucell_->iat2it[na];
            const int iia = ucell_->iat2ia[na];
            const double za = ucell_->atoms[ita].ncpp.zv;
            const double ma = ucell_->atoms[ita].mass;
            const ModuleBase::Vector3<double> &ta = ucell_->atoms[ita].tau[iia];
            for (int nb = 0; nb < nat; ++nb)
            {
                const int itb = ucell_->iat2it[nb];
                const int iib = ucell_->iat2ia[nb];
                const double zb = ucell_->atoms[itb].ncpp.zv;
                const double mb = ucell_->atoms[itb].mass;
                const ModuleBase::Vector3<double> &tb = ucell_->atoms[itb].tau[iib];
                const ModuleBase::Vector3<double> dtau = ta - tb;
                const double pref_pair = ModuleBase::e2 * za * zb / std::sqrt(ma * mb);
                const double pref_diag = ModuleBase::e2 * za * zb / ma;
                H_Ewald_pw::rgen(dtau, rmax_lat0, irr.data(), ucell_->latvec, G, r.data(), r2.data(), mxr, nrm, 0);
                for (int nr = 0; nr < nrm; ++nr)
                {
                    const double rr = std::sqrt(r2[nr]) * ucell_->lat0;
                    const double ar = std::sqrt(ewald_alpha_) * rr;
                    double d2f = 0.0;
                    double df = 0.0;
                    ewald_radial(ar, rr, d2f, df);
                    // recover the integer lattice triplet n of R = r + dtau
                    // (a_i . g_j = delta_ij in ABACUS, so n = R . G^T); the
                    // phase argument stays exact in the small integers
                    const ModuleBase::Vector3<double> R = r[nr] + dtau;
                    const double proj1 = R.x * G.e11 + R.y * G.e12 + R.z * G.e13;
                    const double proj2 = R.x * G.e21 + R.y * G.e22 + R.z * G.e23;
                    const double proj3 = R.x * G.e31 + R.y * G.e32 + R.z * G.e33;
                    const int n1 = static_cast<int>(std::lround(proj1));
                    const int n2 = static_cast<int>(std::lround(proj2));
                    const int n3 = static_cast<int>(std::lround(proj3));
                    const double ph_arg = ModuleBase::TWO_PI * (q_frac.x * n1 + q_frac.y * n2 + q_frac.z * n3);
                    const std::complex<double> phase(std::cos(ph_arg), std::sin(ph_arg));
                    // r in bohr (cartesian) for the tensor contraction
                    const double rx = r[nr].x * ucell_->lat0;
                    const double ry = r[nr].y * ucell_->lat0;
                    const double rz = r[nr].z * ucell_->lat0;
                    const double rv[3] = {rx, ry, rz};
                    for (int da = 0; da < 3; ++da)
                    {
                        for (int db = 0; db < 3; ++db)
                        {
                            const double tensor = d2f * rv[da] * rv[db];
                            dy3(3 * na + da, 3 * nb + db) += pref_pair * phase * tensor;
                            dy3(3 * na + da, 3 * na + db) -= pref_diag * tensor;
                        }
                        const double iso = df;
                        dy3(3 * na + da, 3 * nb + da) += pref_pair * phase * iso;
                        dy3(3 * na + da, 3 * na + da) -= pref_diag * iso;
                    }
                }
            }
        }
    }

    // ---------------- pool reduction ----------------
    // The G-space loops run over the local charge-grid G list and the
    // real-space sum runs only on the rank owning G = 0, so each matrix
    // is rank-local until reduced. The DFPT esolver currently runs
    // single-rank (Autotest.sh convention), but the reduction keeps the
    // ionic term correct under any future multi-rank run and matches the
    // reduce_pool pattern of the electronic accumulation. The guard is
    // compile-time: the serial test targets compile this file without
    // __MPI, and the base-library reduce_pool wrapper is built with MPI,
    // so an unguarded call would hit MPI_Allreduce outside MPI_Init.
#ifdef __MPI
    Parallel_Reduce::reduce_pool(dy1.c, dy1.size);
    Parallel_Reduce::reduce_pool(dy2.c, dy2.size);
    Parallel_Reduce::reduce_pool(dy3.c, dy3.size);
#endif

    // ---------------- assemble ----------------
    // QE d2ionq: dyn = -(dy1 - dy2 + real-space)
    for (int i = 0; i < nat3; ++i)
    {
        for (int j = 0; j < nat3; ++j)
        {
            dyn(i, j) += -(dy1(i, j) - dy2(i, j) + dy3(i, j));
        }
    }
    ModuleBase::timer::end("DFPT_Phon", "ion_ion");
}

} // namespace ModuleDFPT
