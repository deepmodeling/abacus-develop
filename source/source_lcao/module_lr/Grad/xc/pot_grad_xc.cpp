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
            for (int ir = 0;ir < nrxx_;++ir)
            {
                v_eff(0, ir) += ModuleBase::e2 * g.a_s2.at(ir) * rho[0][ir] * rho[0][ir];
            }
        }
        else if (func_type == 2 || func_type == 4)  // GGA or HYB_GGA
        {
            std::vector<ModuleBase::Vector3<double>> drho1(nrxx_);  // transition density gradient
            LR_Util::grad(rho[0], drho1.data(), this->rho_basis_, this->tpiba_);

            std::vector<double> v_tmp(nrxx_, 0.0);

            // 1. the vector under the divergence, accumulated negated so that `grad_dot` yields
            //    $-\nabla\cdot\boldsymbol{E}$. The four $e$ coefficients share the same
            //    $\nabla\rho^{gs}$ direction (see `KernelXC::GxcCoef`), so it is pulled out of
            //    their sum -- exact, and it keeps this bandwidth-bound loop reading 4 doubles per
            //    point instead of 12.
            std::vector<ModuleBase::Vector3<double>> gdot_terms(nrxx_);
            for (int ir = 0;ir < nrxx_;++ir)
            {
                const ModuleBase::Vector3<double>& drho = kxc.drho_gs.at(0).at(ir);  // $\nabla\rho$
                const double s = rho[0][ir];                                // $\rho^1$
                const double t = drho * drho1.at(ir);                       // $\nabla\rho\cdot\nabla\rho^1$
                const double q = drho1.at(ir) * drho1.at(ir);               // $\nabla\rho^1\cdot\nabla\rho^1$
                const double e = g.e_s2.at(ir) * (s * s) + g.e_st.at(ir) * (s * t)
                    + g.e_t2.at(ir) * (t * t) + g.e_q.at(ir) * q;
                gdot_terms[ir] = -(drho * e + drho1.at(ir) * (g.c_s.at(ir) * s + g.c_t.at(ir) * t));
            }
            XC_Functional::grad_dot(gdot_terms.data(), v_tmp.data(), &this->rho_basis_, this->tpiba_);

            // 2. the local terms $A$
            for (int ir = 0;ir < nrxx_;++ir)
            {
                const double s = rho[0][ir];
                const double t = kxc.drho_gs.at(0).at(ir) * drho1.at(ir);
                const double q = drho1.at(ir) * drho1.at(ir);
                v_tmp[ir] += g.a_s2.at(ir) * (s * s) + g.a_st.at(ir) * (s * t)
                    + g.a_t2.at(ir) * (t * t) + g.a_q.at(ir) * q;
            }
            BlasConnector::axpy(nrxx_, ModuleBase::e2, v_tmp.data(), 1, v_eff.c, 1);
        }
        else
        {
            throw std::domain_error("GlobalV::XC_Functional::get_func_type() =" + std::to_string(func_type)
                + " unfinished in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
        }

        ModuleBase::timer::end("PotGradXCLR", "cal_v_eff");
    }

}
