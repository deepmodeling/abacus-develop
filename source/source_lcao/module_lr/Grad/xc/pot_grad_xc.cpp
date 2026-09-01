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

    void PotGradXCLR::Scratch::alloc(const int nrxx, const bool gga)
    {
        // resize() on an already-large vector is a no-op, so only the first call allocates.
        if (static_cast<int>(this->vtmp.size()) < nrxx) { this->vtmp.resize(nrxx); }
        if (gga)
        {
            if (static_cast<int>(this->gdot.size()) < nrxx) { this->gdot.resize(nrxx); }
            if (static_cast<int>(this->drho1.size()) < nrxx) { this->drho1.resize(nrxx); }
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
            scratch().alloc(nrxx_, /*gga=*/true);
            Vec3* const drho1 = scratch().drho1.data();   // transition density gradient
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
}
