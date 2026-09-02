#include "pot_grad_xc.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_lr/potentials/xc_kernel.h"
#include "source_base/timer.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_lcao/module_lr/utils/lr_util.h"
#include "source_lcao/module_lr/utils/lr_util_xc.hpp"
#include <set>
namespace LR
{
    using Vec3 = ModuleBase::Vector3<double>;
    PotGradXCLR::Scratch& PotGradXCLR::scratch()
    {
        static Scratch sc;   // see the comment on `Scratch` in the header
        return sc;
    }

    void PotGradXCLR::Scratch::alloc(const int nrxx, const bool two_channel, const bool gga)
    {
        // resize() on an already-large vector is a no-op, so only the first call allocates.
        if (static_cast<int>(this->vtmp.size()) < nrxx) { this->vtmp.resize(nrxx); }
        if (gga)
        {
            if (static_cast<int>(this->gdot.size()) < nrxx) { this->gdot.resize(nrxx); }
            if (static_cast<int>(this->div.size()) < nrxx) { this->div.resize(nrxx); }
            const int nch = two_channel ? 2 : 1;
            for (int is = 0; is < nch; ++is)
            {
                if (static_cast<int>(this->drho1[is].size()) < nrxx) { this->drho1[is].resize(nrxx); }
            }
        }
    }

    // constructor for exchange-correlation kernel
    PotGradXCLR::PotGradXCLR(const KernelXC& xc_kernel, const ModulePW::PW_Basis& rho_basis, const UnitCell& ucell,
        const int& nrxx, const bool triplet)
        :xc_kernel_components_(xc_kernel), triplet_(triplet),
        PotLRBase(rho_basis, (PARAM.inp.nspin == 1 || (PARAM.inp.nspin == 4 && !PARAM.globalv.domag && !PARAM.globalv.domag_z) ? 1 : 2), nrxx, ucell.tpiba)
    {}

    /// $v^{(2)}(r)=\iint dr'dr''\,g^{xc}(r,r',r'')\rho^1(r')\rho^1(r'')$, i.e. the third functional
    /// derivative of $E_{xc}$ contracted twice with the transition density $\rho^1$ (no factor 1/2).
    ///
    /// All the coefficients and spin sums live in `KernelXC::GxcCoef`, so this is a plain
    /// transcription of the boxed formula and is identical for nspin=1, singlet and triplet.
    ///
    /// Worth stating once: for a GGA, $g^{xc}$ is NOT the whole of $v^{(2)}$. Because
    /// $\sigma=\nabla\rho\cdot\nabla\rho$ is quadratic in the density it has a non-vanishing
    /// *second* derivative along $\rho^1$, which drags the second-order kernels $f^{\rho\sigma}$
    /// and $f^{\sigma\sigma}$ into the answer (the $a_q$, $\boldsymbol{e}_q$, $c_s$, $c_t$ terms).
    void PotGradXCLR::cal_v_eff(double** rho, const UnitCell& ucell, ModuleBase::matrix& v_eff, const std::vector<int>& ispin_op) const
    {
        ModuleBase::TITLE("PotGradXCLR", "cal_v_eff");
        ModuleBase::timer::start("PotGradXCLR", "cal_v_eff");
        const int func_type = XC_Functional::get_func_type();
        const auto& kxc = this->xc_kernel_components_;

        if (kxc.openshell)
        {
            throw std::domain_error("open shell (S2_updown) unfinished in "
                + std::string(__FILE__) + " line " + std::to_string(__LINE__));
        }
        const auto& g = kxc.gxc(this->triplet_);

        if (func_type == 1) // LDA: only the $g^{\rho\rho\rho}$ term survives
        {
            const double* const a_s2 = g.a_s2.data();
            const double* const r1 = rho[0];
            double* const v = v_eff.c;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (int ir = 0;ir < nrxx_;++ir)
            {
                v[ir] += ModuleBase::e2 * a_s2[ir] * r1[ir] * r1[ir];
            }
        }
        else if (func_type == 2 || func_type == 4)  // GGA or HYB_GGA
        {
            scratch().alloc(nrxx_, /*two_channel=*/false, /*gga=*/true);
            Vec3* const drho1 = scratch().drho1[0].data();   // transition density gradient
            LR_Util::grad(rho[0], drho1, this->rho_basis_, this->tpiba_);

            double* const v_tmp = scratch().vtmp.data();
            Vec3* const gdot_terms = scratch().gdot.data();
            const Vec3* const dgs = kxc.drho_gs.at(0).data();
            const double* const r1 = rho[0];
            const double* const e_s2 = g.e_s2.data(); const double* const e_st = g.e_st.data();
            const double* const e_t2 = g.e_t2.data(); const double* const e_q = g.e_q.data();
            const double* const c_s = g.c_s.data();   const double* const c_t = g.c_t.data();
            const double* const a_s2 = g.a_s2.data(); const double* const a_st = g.a_st.data();
            const double* const a_t2 = g.a_t2.data(); const double* const a_q = g.a_q.data();

            // 1. the vector under the divergence, accumulated negated so that `grad_dot` yields
            //    $-\nabla\cdot\boldsymbol{E}$. The four $e$ coefficients share the same
            //    $\nabla\rho^{gs}$ direction (see `KernelXC::GxcCoef`), so it is pulled out of
            //    their sum -- exact, and it keeps this bandwidth-bound loop reading 4 doubles per
            //    point instead of 12.
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (int ir = 0;ir < nrxx_;++ir)
            {
                const Vec3& drho = dgs[ir];                     // $\nabla\rho$
                const double s = r1[ir];                        // $\rho^1$
                const double t = drho * drho1[ir];              // $\nabla\rho\cdot\nabla\rho^1$
                const double q = drho1[ir] * drho1[ir];         // $\nabla\rho^1\cdot\nabla\rho^1$
                const double e = e_s2[ir] * (s * s) + e_st[ir] * (s * t)
                    + e_t2[ir] * (t * t) + e_q[ir] * q;
                gdot_terms[ir] = -(drho * e + drho1[ir] * (c_s[ir] * s + c_t[ir] * t));
            }
            XC_Functional::grad_dot(gdot_terms, v_tmp, &this->rho_basis_, this->tpiba_);

            // 2. the local terms $A$
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (int ir = 0;ir < nrxx_;++ir)
            {
                const double s = r1[ir];
                const double t = dgs[ir] * drho1[ir];
                const double q = drho1[ir] * drho1[ir];
                v_tmp[ir] += a_s2[ir] * (s * s) + a_st[ir] * (s * t)
                    + a_t2[ir] * (t * t) + a_q[ir] * q;
            }
            BlasConnector::axpy(nrxx_, ModuleBase::e2, v_tmp, 1, v_eff.c, 1);
        }
        else
        {
            throw std::domain_error("GlobalV::XC_Functional::get_func_type() =" + std::to_string(func_type)
                + " unfinished in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
        }

        ModuleBase::timer::end("PotGradXCLR", "cal_v_eff");
    }


    /// $v^{(2)}_\tau=A_\tau-\nabla\cdot\boldsymbol{E}_\tau$ for the open-shell case, contracted
    /// straight out of the raw libxc arrays (there is no useful pre-contraction: the free spin
    /// $\tau$ stays open, so a `GxcCoef`-style cache would cost ~117 doubles per grid point
    /// against the 35 the third-order arrays already occupy).
    ///
    /// With $s_\sigma=\rho^1_\sigma$, $t_{ab}=\nabla\rho_a\cdot\nabla\rho^1_b$ (NOT symmetric)
    /// and $q_{ab}=\nabla\rho^1_a\cdot\nabla\rho^1_b$, the $\lambda$-derivatives of libxc's
    /// three sigma variables are
    ///     $S=(2t_{uu},\;t_{ud}+t_{du},\;2t_{dd})$,  $Q=(2q_{uu},\;2q_{ud},\;2q_{dd})$,
    /// and with $D=\sum_\sigma s_\sigma\partial_{\rho_\sigma}+\sum_a S_a\partial_{\sigma_a}$,
    ///     $A_\tau = D^2 e^{\rho_\tau} + \sum_a Q_a e^{\rho_\tau\sigma_a}$,
    ///     $u'_a  = D\,e^{\sigma_a}$,  $u''_a = D^2 e^{\sigma_a}+\sum_b Q_b e^{\sigma_a\sigma_b}$,
    ///     $\boldsymbol{E}_\tau=\sum_a\theta^\tau_a
    ///        \big(u''_a\,\nabla\rho_{c(\tau,a)} + 2u'_a\,\nabla\rho^1_{c(\tau,a)}\big)$.
    /// The channel selector $c(\tau,a)$ is what produces the closed-shell $\tilde\theta$: for the
    /// triplet $\nabla\rho^1_d=-\nabla\rho^1_u$, which is invisible in any singlet-only test.
    void PotGradXCLR::cal_v_eff_openshell(const double* const* const rho1, const UnitCell& ucell,
        ModuleBase::matrix& v_eff, const int tau) const
    {
        ModuleBase::TITLE("PotGradXCLR", "cal_v_eff_openshell");
        ModuleBase::timer::start("PotGradXCLR", "cal_v_eff_openshell");
        using namespace LR::libxc_idx;
        const int func_type = XC_Functional::get_func_type();
        const auto& kxc = this->xc_kernel_components_;
        assert(tau == 0 || tau == 1);
        if (func_type != 1 && func_type != 2 && func_type != 4)
        {
            throw std::domain_error("PotGradXCLR: func_type = " + std::to_string(func_type)
                + " (meta-GGA) is not supported, in " + std::string(__FILE__));
        }
        const std::vector<double>& v2rs = kxc.v2rhosigma;
        const std::vector<double>& v2s2 = kxc.v2sigma2;
        const std::vector<double>& v3r3 = kxc.v3rho3;
        const std::vector<double>& v3r2s = kxc.v3rho2sigma;
        const std::vector<double>& v3rs2 = kxc.v3rhosigma2;
        const std::vector<double>& v3s3 = kxc.v3sigma3;

        if (func_type == 1)   // LDA: only $g^{\rho\rho\rho}$ survives
        {
            const double* const r1u = rho1[0]; const double* const r1d = rho1[1];
            const double* const g3 = v3r3.data();
            double* const v = v_eff.c;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (int ir = 0;ir < nrxx_;++ir)
            {
                const double s[2] = { r1u[ir], r1d[ir] };
                double a = 0.;
                for (int s0 = 0;s0 < 2;++s0) {
                    for (int s1 = 0;s1 < 2;++s1) { a += s[s0] * s[s1] * g3[ir * 4 + r3(tau, s0, s1)]; } }
                v[ir] += ModuleBase::e2 * a;
            }
            ModuleBase::timer::end("PotGradXCLR", "cal_v_eff_openshell");
            return;
        }

        // GGA / HYB_GGA
        scratch().alloc(nrxx_, /*two_channel=*/true, /*gga=*/true);
        Vec3* const drho1[2] = { scratch().drho1[0].data(), scratch().drho1[1].data() };
        for (int is : {0, 1}) { LR_Util::grad(rho1[is], drho1[is], this->rho_basis_, this->tpiba_); }

        double* const v_tmp = scratch().vtmp.data();
        Vec3* const gdot_terms = scratch().gdot.data();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ir = 0;ir < nrxx_;++ir)
        {
            const int o4 = ir * 4, o6 = ir * 6, o9 = ir * 9, o10 = ir * 10, o12 = ir * 12;
            const ModuleBase::Vector3<double> drho[2] = { kxc.drho_gs[0][ir], kxc.drho_gs[1][ir] };
            const ModuleBase::Vector3<double> dr1[2] = { drho1[0][ir], drho1[1][ir] };
            const double s[2] = { rho1[0][ir], rho1[1][ir] };

            // $t_{ab}=\nabla\rho_a\cdot\nabla\rho^1_b$, then $S_a=\mathrm{d}\sigma_a/\mathrm{d}\lambda$
            const double t00 = drho[0] * dr1[0], t01 = drho[0] * dr1[1];
            const double t10 = drho[1] * dr1[0], t11 = drho[1] * dr1[1];
            const double S[3] = { 2. * t00, t01 + t10, 2. * t11 };
            // $Q_a=\mathrm{d}^2\sigma_a/\mathrm{d}\lambda^2$
            const double Q[3] = { 2. * (dr1[0] * dr1[0]), 2. * (dr1[0] * dr1[1]), 2. * (dr1[1] * dr1[1]) };

            // ---- the local part $A_\tau$ ----
            double A = 0.;
            for (int s0 = 0;s0 < 2;++s0) {
                for (int s1 = 0;s1 < 2;++s1) { A += s[s0] * s[s1] * v3r3[o4 + r3(tau, s0, s1)]; } }
            for (int s0 = 0;s0 < 2;++s0) {
                for (int a = 0;a < 3;++a) { A += 2. * s[s0] * S[a] * v3r2s[o9 + r2s(tau, s0, a)]; } }
            for (int a = 0;a < 3;++a) {
                for (int b = 0;b < 3;++b) { A += S[a] * S[b] * v3rs2[o12 + rs2(tau, a, b)]; } }
            for (int a = 0;a < 3;++a) { A += Q[a] * v2rs[o6 + rs(tau, a)]; }
            v_tmp[ir] = A;

            // ---- the divergence part $\boldsymbol{E}_\tau$ ----
            // $u'_a$ and $u''_a$ do not depend on $\tau$; only the $\theta$ weights and the
            // channel selector below do.
            ModuleBase::Vector3<double> E(0., 0., 0.);
            for (int a = 0;a < 3;++a)
            {
                const double th = theta[tau][a];
                if (th == 0.) { continue; }
                const int c = chan[tau][a];
                double up = 0., upp = 0.;
                for (int s0 = 0;s0 < 2;++s0) { up += s[s0] * v2rs[o6 + rs(s0, a)]; }
                for (int b = 0;b < 3;++b) { up += S[b] * v2s2[o6 + p2[a][b]]; }

                for (int s0 = 0;s0 < 2;++s0) {
                    for (int s1 = 0;s1 < 2;++s1) { upp += s[s0] * s[s1] * v3r2s[o9 + r2s(s0, s1, a)]; } }
                for (int s0 = 0;s0 < 2;++s0) {
                    for (int b = 0;b < 3;++b) { upp += 2. * s[s0] * S[b] * v3rs2[o12 + rs2(s0, a, b)]; } }
                for (int b = 0;b < 3;++b) {
                    for (int cc = 0;cc < 3;++cc) { upp += S[b] * S[cc] * v3s3[o10 + p3[a][b][cc]]; } }
                for (int b = 0;b < 3;++b) { upp += Q[b] * v2s2[o6 + p2[a][b]]; }

                E += th * (drho[c] * upp + dr1[c] * (2. * up));
            }
            gdot_terms[ir] = -E;   // `grad_dot` then yields $-\nabla\cdot\boldsymbol{E}_\tau$
        }
        double* const div = scratch().div.data();
        XC_Functional::grad_dot(gdot_terms, div, &this->rho_basis_, this->tpiba_);
        double* const vout = v_eff.c;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ir = 0;ir < nrxx_;++ir) { vout[ir] += ModuleBase::e2 * (v_tmp[ir] + div[ir]); }
        ModuleBase::timer::end("PotGradXCLR", "cal_v_eff_openshell");
    }
}
