// Ewald ion-ion part of DFPT_Phon::ion_ion, split out of
// dfpt_phon.cpp: the screening-alpha search, the reciprocal-space
// Poisson sums and the real-space erfc force constants. All formulas
// are moved verbatim from the original ion_ion body.

#include "dfpt_phon.h"

#include "source_base/constants.h"
#include "source_base/tool_quit.h"
#include "source_base/truncated_func.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/unitcell.h"

#include <cmath>
#include <complex>
#include <vector>

#include "source_base/timer.h"
#include "source_base/tool_title.h"

namespace ModuleDFPT
{

namespace
{

/// converged screening alpha of the Ewald split, searched so that the
/// G-sum tail is bounded inside the rho grid (returns the alpha; the
/// caller stores it and derives the real-space cutoff)
double ewald_alpha_search(double charge, double ggecut, double tpiba2)
{
    const double alpha_init = 1.1;        ///< empirical parameter: initial Ewald screening-alpha guess
    const double alpha_shrink = 0.9;      ///< empirical parameter: per-iteration alpha shrink factor
    const double alpha_min = 1.0e-4;      ///< empirical parameter: alpha floor of the search
    const double ewald_tail_tol = 1.0e-6; ///< empirical parameter: accepted G-sum tail bound
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
        upperbound = 2.0 * charge * charge * std::sqrt(2.0 * alpha / ModuleBase::TWO_PI)
                     * ModuleBase::truncated_erfc(std::sqrt(ggecut * tpiba2 / 4.0 / alpha));
    } while (upperbound > ewald_tail_tol);
    return alpha;
}

/// sq/s0 kernels of the self-image phase difference: sq sums
/// (G+q)(G+q)/|G+q|^2 exp(-|G+q|^2/4a) over all grid G (the G = 0 member
/// contributes through w = q) and s0 the same kernel at q = 0
void g_self_accum(const UnitCell& ucell,
                  const ModulePW::PW_Basis& pw_rho,
                  const ModuleBase::Vector3<double>& q_cart,
                  double alpha,
                  double (&sq)[3][3],
                  double (&s0)[3][3])
{
    const double w2_floor = 1.0e-12; ///< empirical parameter: |G+q|^2 zero-shell guard (1/lat0^2)
    const double g2_floor = 1.0e-12; ///< empirical parameter: |G|^2 zero-shell guard (1/lat0^2)
    for (int ig = 0; ig < pw_rho.npw; ++ig)
    {
        const ModuleBase::Vector3<double>& gcart = pw_rho.gcar[ig];
        const ModuleBase::Vector3<double> w = gcart + q_cart;
        const double w2 = w * w;
        const double g2 = gcart * gcart;
        if (w2 < w2_floor)
        {
            // G + q = 0 (only possible at q = 0 with G = 0): excluded, as in
            // the q = 0 G part below; its isotropic delta/3 limit belongs to
            // the direction-averaged q -> 0 behavior, not the exact q = 0
            // matrix
            continue;
        }
        const double w2_bohr = w2 * ucell.tpiba2;
        const double gauss = ModuleBase::truncated_exp(-w2_bohr / (4.0 * alpha));
        for (int da = 0; da < 3; ++da)
        {
            for (int db = 0; db < 3; ++db)
            {
                sq[da][db] += w[da] * w[db] / w2 * gauss;
            }
        }
        if (g2 > g2_floor)
        {
            const double gauss_g = ModuleBase::truncated_exp(-g2 * ucell.tpiba2 / (4.0 * alpha));
            for (int da = 0; da < 3; ++da)
            {
                for (int db = 0; db < 3; ++db)
                {
                    s0[da][db] += gcart[da] * gcart[db] / g2 * gauss_g;
                }
            }
        }
    }
}

/// reciprocal-space Poisson pair sum over atom pairs (ia != ib) with the
/// phase-free Gamma on-site diagonal of the same-atom block
void g_pair_accum(const UnitCell& ucell,
                  const ModulePW::PW_Basis& pw_rho,
                  const ModuleBase::Vector3<double>& q_cart,
                  double alpha,
                  ModuleBase::ComplexMatrix& dyn)
{
    const double w2_floor = 1.0e-12; ///< empirical parameter: |G+q|^2 zero-shell guard (1/lat0^2)
    const double g2_floor = 1.0e-12; ///< empirical parameter: |G|^2 zero-shell guard (1/lat0^2)
    const int nat = ucell.nat;
    for (int ig = 0; ig < pw_rho.npw; ++ig)
    {
        const ModuleBase::Vector3<double>& gcart = pw_rho.gcar[ig];
        const ModuleBase::Vector3<double> w = gcart + q_cart;
        const double w2 = w * w;
        const double g2 = gcart * gcart;
        if (w2 < w2_floor)
        {
            // G + q = 0 (only possible at q = 0 with G = 0): excluded
            continue;
        }
        const double w2_bohr = w2 * ucell.tpiba2;
        const double gauss = ModuleBase::truncated_exp(-w2_bohr / (4.0 * alpha));
        double gauss_g = 0.0;
        if (g2 > g2_floor)
        {
            gauss_g = ModuleBase::truncated_exp(-g2 * ucell.tpiba2 / (4.0 * alpha));
        }
        for (int ia = 0; ia < nat; ++ia)
        {
            const int ita = ucell.iat2it[ia];
            const int iia = ucell.iat2ia[ia];
            const double za = ucell.atoms[ita].ncpp.zv;
            const double ma = ucell.atoms[ita].mass;
            const ModuleBase::Vector3<double>& ta = ucell.atoms[ita].tau[iia];
            for (int ib = 0; ib < nat; ++ib)
            {
                if (ib == ia)
                {
                    continue;
                }
                const int itb = ucell.iat2it[ib];
                const int iib = ucell.iat2ia[ib];
                const double zb = ucell.atoms[itb].ncpp.zv;
                const double mb = ucell.atoms[itb].mass;
                const ModuleBase::Vector3<double>& tb = ucell.atoms[itb].tau[iib];
                const double arg = ModuleBase::TWO_PI * (w * (ta - tb));
                const std::complex<double> phase(std::cos(arg), std::sin(arg));
                const double pref = ModuleBase::FOUR_PI / ucell.omega * za * zb * ModuleBase::e2 * gauss
                                    / (std::sqrt(ma * mb) * w2);
                // Gamma-phase on-site piece (G-only kernel, G != 0)
                std::complex<double> phase0(1.0, 0.0);
                double pref0 = 0.0;
                if (g2 > g2_floor)
                {
                    const double arg0 = ModuleBase::TWO_PI * (gcart * (ta - tb));
                    phase0 = std::complex<double>(std::cos(arg0), std::sin(arg0));
                    pref0 = ModuleBase::FOUR_PI / ucell.omega * za * zb * ModuleBase::e2 * gauss_g
                            / (std::sqrt(ma * mb) * g2);
                }
                for (int da = 0; da < 3; ++da)
                {
                    for (int db = 0; db < 3; ++db)
                    {
                        const std::complex<double> elem = pref * w[da] * w[db] * phase;
                        dyn(3 * ia + da, 3 * ib + db) += elem;
                        // on-site diagonal: phase-free (Gamma) accumulation,
                        // Phi_ii = -Phi_ij => -sqrt(Mb/Ma) on the pair term
                        dyn(3 * ia + da, 3 * ia + db) -= pref0 * gcart[da] * gcart[db] * phase0 * std::sqrt(mb / ma);
                    }
                }
            }
        }
    }
}

/// d^2/dR_a dR_b [ erfc(sqrt(alpha) R) / R ] tensor (bohr^-3), shared by
/// the self-image and the pair real-space sums; validated against central
/// finite differences of the erfc-split Ewald energy
void ewald_h_ab(const ModuleBase::Vector3<double>& r, double alpha, double (&h)[3][3])
{
    const double r2 = r * r;
    const double rlen = std::sqrt(r2);
    const double sar = std::sqrt(alpha);
    const double e2a = ModuleBase::truncated_exp(-alpha * r2);
    const double f = 2.0 * sar / std::sqrt(ModuleBase::PI) * e2a;
    const double er = ModuleBase::truncated_erfc(sar * rlen);
    for (int da = 0; da < 3; ++da)
    {
        for (int db = 0; db < 3; ++db)
        {
            const double delta = (da == db) ? 1.0 : 0.0;
            h[da][db] = er * (3.0 * r[da] * r[db] - delta * r2) / (rlen * r2 * r2)
                        + f
                              * (2.0 * alpha * r[da] * r[db] / r2
                                 + 3.0 * r[da] * r[db] / (r2 * r2) - delta / r2);
        }
    }
}

/// ranges of the lattice-vector shells (rows of latvec are the lattice
/// translations in lat0 units) covering the real-space cutoff
void real_shell_ranges(const ModuleBase::Matrix3& latvec, double lat0, double rcut, int (&nmax)[3])
{
    const double row_e[3][3] = {{latvec.e11, latvec.e12, latvec.e13},
                                {latvec.e21, latvec.e22, latvec.e23},
                                {latvec.e31, latvec.e32, latvec.e33}};
    for (int d = 0; d < 3; ++d)
    {
        const ModuleBase::Vector3<double> a1(row_e[d][0], row_e[d][1], row_e[d][2]);
        const double len = std::sqrt(a1 * a1) * lat0; // bohr
        nmax[d] = static_cast<int>(std::ceil(rcut / len)) + 1;
    }
}

/// real-space self-image phase difference of the on-site i-i block: the
/// on-site i-i energy is L-independent while the cross-cell i-i force
/// constants carry e^{i2pi q.L}, so D_ii receives
///   -(Za^2 e2/Ma) sum_{L!=0} h_erfc(L) (e^{i2pi q.L} - 1);
/// the imaginary part cancels over the +-L symmetric sphere (h is even)
/// and L = 0 carries e^{i2pi q.0} - 1 = 0
void real_self_images(const UnitCell& ucell,
                      const ModuleBase::Vector3<double>& q_frac,
                      double lat0,
                      double alpha,
                      double rcut,
                      const int (&nmax)[3],
                      ModuleBase::ComplexMatrix& dyn)
{
    const ModuleBase::Matrix3& latvec = ucell.latvec;
    for (int ia = 0; ia < ucell.nat; ++ia)
    {
        const int ita = ucell.iat2it[ia];
        const int iia = ucell.iat2ia[ia];
        const double za = ucell.atoms[ita].ncpp.zv;
        const double ma = ucell.atoms[ita].mass;
        const double f2 = za * za * ModuleBase::e2 / ma;
        for (int n1 = -nmax[0]; n1 <= nmax[0]; ++n1)
        {
            for (int n2 = -nmax[1]; n2 <= nmax[1]; ++n2)
            {
                for (int n3 = -nmax[2]; n3 <= nmax[2]; ++n3)
                {
                    if (n1 == 0 && n2 == 0 && n3 == 0)
                    {
                        continue;
                    }
                    const ModuleBase::Vector3<double> lvec(n1 * latvec.e11 + n2 * latvec.e21 + n3 * latvec.e31,
                                                           n1 * latvec.e12 + n2 * latvec.e22 + n3 * latvec.e32,
                                                           n1 * latvec.e13 + n2 * latvec.e23 + n3 * latvec.e33);
                    const ModuleBase::Vector3<double> r = lvec * lat0;
                    const double r2 = r * r;
                    if (r2 > rcut * rcut)
                    {
                        continue;
                    }
                    const double ph_arg = ModuleBase::TWO_PI * (q_frac.x * n1 + q_frac.y * n2 + q_frac.z * n3);
                    const double wcos = std::cos(ph_arg) - 1.0;
                    double h[3][3];
                    ewald_h_ab(r, alpha, h);
                    for (int da = 0; da < 3; ++da)
                    {
                        for (int db = 0; db < 3; ++db)
                        {
                            dyn(3 * ia + da, 3 * ia + db) -= f2 * h[da][db] * wcos;
                        }
                    }
                }
            }
        }
    }
}

/// real-space pair sum over lattice shells for ia != ib with the phase-free
/// on-site diagonal Phi_ii^R = -sqrt(Mb/Ma) times the pair term
void real_pair_sum(const UnitCell& ucell,
                   const ModuleBase::Vector3<double>& q_frac,
                   double lat0,
                   double alpha,
                   double rcut,
                   const int (&nmax)[3],
                   ModuleBase::ComplexMatrix& dyn)
{
    const ModuleBase::Matrix3& latvec = ucell.latvec;
    for (int ia = 0; ia < ucell.nat; ++ia)
    {
        const int ita = ucell.iat2it[ia];
        const int iia = ucell.iat2ia[ia];
        const double za = ucell.atoms[ita].ncpp.zv;
        const double ma = ucell.atoms[ita].mass;
        for (int ib = 0; ib < ucell.nat; ++ib)
        {
            if (ib == ia)
            {
                continue;
            }
            const int itb = ucell.iat2it[ib];
            const int iib = ucell.iat2ia[ib];
            const double zb = ucell.atoms[itb].ncpp.zv;
            const double mb = ucell.atoms[itb].mass;
            const ModuleBase::Vector3<double> dt = ucell.atoms[itb].tau[iib] - ucell.atoms[ita].tau[iia];
            const double zab2 = za * zb * ModuleBase::e2 / std::sqrt(ma * mb);
            for (int n1 = -nmax[0]; n1 <= nmax[0]; ++n1)
            {
                for (int n2 = -nmax[1]; n2 <= nmax[1]; ++n2)
                {
                    for (int n3 = -nmax[2]; n3 <= nmax[2]; ++n3)
                    {
                        const ModuleBase::Vector3<double> lvec(n1 * latvec.e11 + n2 * latvec.e21 + n3 * latvec.e31,
                                                               n1 * latvec.e12 + n2 * latvec.e22 + n3 * latvec.e32,
                                                               n1 * latvec.e13 + n2 * latvec.e23 + n3 * latvec.e33);
                        ModuleBase::Vector3<double> r = (lvec + dt) * lat0; // bohr
                        const double r2 = r * r;
                        if (r2 > rcut * rcut)
                        {
                            continue;
                        }
                        const double ph_arg = ModuleBase::TWO_PI * (q_frac.x * n1 + q_frac.y * n2 + q_frac.z * n3);
                        const std::complex<double> phase(std::cos(ph_arg), std::sin(ph_arg));
                        double h[3][3];
                        ewald_h_ab(r, alpha, h);
                        for (int da = 0; da < 3; ++da)
                        {
                            for (int db = 0; db < 3; ++db)
                            {
                                dyn(3 * ia + da, 3 * ib + db) -= zab2 * h[da][db] * phase;
                                // on-site diagonal Phi_ii^R = sum_{j != i}
                                // Z_iZ_j sum_L h(r_ij + L): phase-free (both
                                // derivatives act on tau_a in cell 0), i.e.
                                // -sqrt(Mb/Ma) times the pair term
                                dyn(3 * ia + da, 3 * ia + db) += zab2 * std::sqrt(mb / ma) * h[da][db];
                            }
                        }
                    }
                }
            }
        }
    }
}

} // namespace

void DFPT_Phon::ion_ion(const ModuleBase::Vector3<double>& q_frac, ModuleBase::ComplexMatrix& dyn)
{
    ModuleBase::TITLE("DFPT_Phon", "ion_ion");
    ModuleBase::timer::start("DFPT_Phon", "ion_ion");
    const double lat0 = ucell_->lat0;

    // total ionic charge
    double charge = 0.0;
    for (int it = 0; it < ucell_->ntype; ++it)
    {
        charge += ucell_->atoms[it].na * ucell_->atoms[it].ncpp.zv;
    }

    // choose the screening alpha so that the G-sum tail is converged inside
    // the rho grid (the erfc envelope bounds the exp(-G^2/4alpha) tail)
    ewald_alpha_ = ewald_alpha_search(charge, pw_rho_->ggecut, ucell_->tpiba2);
    const double ewald_rcut_factor = 6.0; ///< empirical parameter: real-space cutoff in 1/sqrt(alpha)
    // erfc(alpha R) < 1e-16 well inside 6/sqrt(alpha)
    ewald_rcut_ = ewald_rcut_factor / std::sqrt(ewald_alpha_);

    const ModuleBase::Vector3<double> q_cart = q_frac * ucell_->G;

    // ---------------- reciprocal-space part ----------------
    // Poisson pair identity (validated against direct sums):
    //   sum_L h(R) e^{i2pi q.L} = sum_L h_erfc(R) e^{i2pi q.L}
    //     + (4pi/Omega) sum_{|G+q|>0} (G+q)_a (G+q)_b / |G+q|^2
    //       exp(-|G+q|^2/4a) e^{i2pi (G+q).(tau_a-tau_b)}
    // so the G part enters D with the + sign while the erfc part carries -.
    // The on-site diagonal (both second derivatives act on tau_a in cell 0)
    // is phase-free: it is accumulated from Gamma-phase (G-only) pair terms
    // as -sqrt(Mb/Ma) times the pair element. sq/s0 accumulate the self-image
    // phase difference of the same-atom images (validated element-wise
    // against finite differences of the erfc-split Ewald energy in a
    // q-commensurate supercell):
    //   D_ii(q) - D_ii(0) = (Za^2 e2 / Ma) [ sum_{L!=0} h(L)(1 - cos(2pi q.L))
    //     + (4pi/Omega)(sq - s0) ],
    // where sq/s0 are the kernels collected by g_self_accum. The alpha
    // independence of this combination was verified numerically; at q = 0
    // both differences vanish and the acoustic sum rule holds exactly by
    // construction.
    double sq[3][3] = {{0.0}};
    double s0[3][3] = {{0.0}};
    g_self_accum(*ucell_, *pw_rho_, q_cart, ewald_alpha_, sq, s0);
    g_pair_accum(*ucell_, *pw_rho_, q_cart, ewald_alpha_, dyn);
    // self-image G-space phase difference on the diagonal
    for (int ia = 0; ia < ucell_->nat; ++ia)
    {
        const int ita = ucell_->iat2it[ia];
        const double za = ucell_->atoms[ita].ncpp.zv;
        const double ma = ucell_->atoms[ita].mass;
        const double f2 = za * za * ModuleBase::e2 / ma;
        for (int da = 0; da < 3; ++da)
        {
            for (int db = 0; db < 3; ++db)
            {
                dyn(3 * ia + da, 3 * ia + db) += f2 * ModuleBase::FOUR_PI / ucell_->omega * (sq[da][db] - s0[da][db]);
            }
        }
    }

    // ---------------- real-space part ----------------
    // h_ab(R) = d^2/dR_a dR_b [ erfc(sqrt(alpha) R) / R ]
    //         = erfc(sqrt(alpha) R) (3 Ra Rb - delta R^2)/R^5
    //         + (2 sqrt(alpha)/sqrt(pi)) e^{-alpha R^2}
    //           [ 2 alpha Ra Rb/R^2 + 3 Ra Rb/R^4 - delta/R^2 ]
    // D^R_ab = -(1/sqrt(MaMb)) ZaZb e2 h(R = tau_b + l - tau_a) e^{i2pi q.l}
    int nmax[3] = {0, 0, 0};
    real_shell_ranges(ucell_->latvec, lat0, ewald_rcut_, nmax);
    real_self_images(*ucell_, q_frac, lat0, ewald_alpha_, ewald_rcut_, nmax, dyn);
    real_pair_sum(*ucell_, q_frac, lat0, ewald_alpha_, ewald_rcut_, nmax, dyn);

    // The Gaussian self constant -Z^2 sqrt(2 alpha/pi) and the h_erf contact
    // -4 alpha^{3/2}/(3 sqrt(pi)) delta_ab are tau-independent and cancel in
    // the (e^{i2pi q.L} - 1) differences; the diagonal is carried by the
    // phase-free cross-atom accumulation plus the self-image phase terms
    // (both G and R pieces above). At q = 0 all phase differences vanish and
    // the acoustic sum rule holds exactly by construction.
    ModuleBase::timer::end("DFPT_Phon", "ion_ion");
}

} // namespace ModuleDFPT
