#include "pot_hxc_lrtd.h"
#include "module_elecstate/potentials/pot_base.h"
#include "module_elecstate/potentials/H_Hartree_pw.h"
#include "module_base/timer.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include <set>
#include "module_lr/utils/lr_util.h"
namespace LR
{
    // constructor for exchange-correlation kernel
    PotHxcLR::PotHxcLR(const std::string& xc_kernel_in, const ModulePW::PW_Basis* rho_basis_in, const UnitCell* ucell_in, const Charge* chg_gs/*ground state*/)
        :xc_kernel(xc_kernel_in)
    {
        std::cout << "xc_kernel_in: " << xc_kernel_in << std::endl;
        this->rho_basis_ = rho_basis_in;
        // this->dynamic_mode = true;
        // this->fixed_mode = false;
        this->nrxx = chg_gs->nrxx;
        this->nspin = (GlobalV::NSPIN == 1 || (GlobalV::NSPIN == 4 && !GlobalV::DOMAG && !GlobalV::DOMAG_Z)) ? 1 : 2;

        //one-time init and cal kernels
        // this->xc_kernel_components_.resize(1, nullptr);
        // this->xc_kernel_components_[0] = new KernelHartree(rho_basis_in);
        this->pot_hartree = new elecstate::PotHartree(this->rho_basis_);
        std::set<std::string> local_xc = { "lda", "pbe", "hse" };
        if (local_xc.find(this->xc_kernel) != local_xc.end())
        {
            XC_Functional::set_xc_type(this->xc_kernel);    // for hse, (1-alpha) and omega are set here
            this->xc_kernel_components_.cal_kernel(chg_gs, ucell_in, this->nspin);
        }
    }

    void PotHxcLR::cal_v_eff(double** rho, const UnitCell* ucell, ModuleBase::matrix& v_eff)
    {
        ModuleBase::TITLE("PotHxcLR", "cal_v_eff");
        ModuleBase::timer::tick("PotHxcLR", "cal_v_eff");
        auto& fxc = this->xc_kernel_components_;
#ifdef USE_LIBXC
        const int nspin = v_eff.nr;
        v_eff += elecstate::H_Hartree_pw::v_hartree(*ucell, const_cast<ModulePW::PW_Basis*>(this->rho_basis_), v_eff.nr, rho);
        if (xc_kernel == "rpa" || xc_kernel == "hf") { return; }
        else if (XC_Functional::get_func_type() == 1)//LDA
        {
            if (1 == nspin)// for LDA-spin0, just fxc*rho where fxc=v2rho2; for GGA, v2rho2 has been replaced by the true fxc
            {
                for (int ir = 0;ir < nrxx;++ir) { v_eff(0, ir) += ModuleBase::e2 * fxc.get_kernel("v2rho2").at(ir) * rho[0][ir]; }
            }
            else if (2 == nspin)
            {
                for (int ir = 0;ir < nrxx;++ir)
                {
                    const int irs0 = 2 * ir;
                    const int irs1 = irs0 + 1;
                    const int irs2 = irs0 + 2;
                    v_eff(0, ir) += ModuleBase::e2 * fxc.get_kernel("v2rho2").at(irs0) * rho[0][ir]
                        + fxc.get_kernel("v2rho2").at(irs1) * rho[1][ir];
                    v_eff(1, ir) += ModuleBase::e2 * fxc.get_kernel("v2rho2").at(irs1) * rho[0][ir]
                        + fxc.get_kernel("v2rho2").at(irs2) * rho[1][ir];
                }
            }
            else  //remain for spin 4
            {
                throw std::domain_error("nspin =" + std::to_string(nspin)
                    + " unfinished in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
            }
        }

        else if (XC_Functional::get_func_type() == 2 || XC_Functional::get_func_type() == 4)    // GGA or HYB_GGA
        {
            if (1 == nspin)
            {
                std::vector<ModuleBase::Vector3<double>> drho(nrxx);    // transition density gradient
                LR_Util::grad(rho[0], drho.data(), *(this->rho_basis_), ucell->tpiba);
                // test: output drho
                double thr = 1e-1;
                auto out_thr = [this, &thr](const double* v) {
                    for (int ir = 0;ir < nrxx;++ir) if (std::abs(v[ir]) > thr) std::cout << v[ir] << " ";
                    std::cout << std::endl;};
                auto out_thr3 = [this, &thr](const std::vector<ModuleBase::Vector3<double>>& v) {
                    for (int ir = 0;ir < nrxx;++ir) if (std::abs(v.at(ir).x) > thr) std::cout << v.at(ir).x << " ";
                    std::cout << std::endl;
                    for (int ir = 0;ir < nrxx;++ir) if (std::abs(v.at(ir).y) > thr) std::cout << v.at(ir).y << " ";
                    std::cout << std::endl;
                    for (int ir = 0;ir < nrxx;++ir) if (std::abs(v.at(ir).z) > thr) std::cout << v.at(ir).z << " ";
                    std::cout << std::endl;};

                std::vector<double> vxc_tmp(nrxx, 0.0);

                //1. $\partial E/\partial\rho = 2f^{\rho\sigma}*\nabla\rho*\rho_1+4f^{\sigma\sigma}\nabla\rho(\nabla\rho\cdot\nabla\rho_1)+2v^\sigma\nabla\rho_1$
                std::vector<ModuleBase::Vector3<double>> e_drho(nrxx);
                for (int ir = 0;ir < nrxx;++ir)
                {
                    e_drho[ir] = -(fxc.get_grad_kernel("2_v2rhosigma_drho").at(ir) * rho[0][ir]
                        + fxc.get_grad_kernel("4_v2sigma2_drho").at(ir) * (fxc.get_grad_kernel("drho_gs").at(ir) * drho.at(ir))
                        + drho.at(ir) * fxc.get_kernel("vsigma").at(ir) * 2.);
                }
                XC_Functional::grad_dot(e_drho.data(), vxc_tmp.data(), this->rho_basis_, ucell->tpiba);

                // 2. $f^{\rho\rho}\rho_1+2f^{\rho\sigma}\nabla\rho\cdot\nabla\rho_1$
                for (int ir = 0;ir < nrxx;++ir)
                {
                    vxc_tmp[ir] += (fxc.get_kernel("v2rho2").at(ir) * rho[0][ir]
                        + fxc.get_grad_kernel("2_v2rhosigma_drho").at(ir) * drho.at(ir));
                }
                BlasConnector::axpy(nrxx, ModuleBase::e2, vxc_tmp.data(), 1, v_eff.c, 1);
            }
            else if (2 == nspin)    // wrong for GGA, to be fixed
            {
                for (int ir = 0;ir < nrxx;++ir)
                {
                    const int irs0 = 2 * ir;
                    const int irs1 = irs0 + 1;
                    const int irs2 = irs0 + 2;
                    v_eff(0, ir) += ModuleBase::e2 * fxc.get_kernel("v2rho2").at(irs0) * rho[0][ir]
                        + fxc.get_kernel("v2rho2").at(irs1) * rho[1][ir];
                    v_eff(1, ir) += ModuleBase::e2 * fxc.get_kernel("v2rho2").at(irs1) * rho[0][ir]
                        + fxc.get_kernel("v2rho2").at(irs2) * rho[1][ir];
                }
            }
            else  //remain for spin 4
            {
                throw std::domain_error("nspin =" + std::to_string(nspin)
                    + " unfinished in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
            }
        }
        else
#endif
        {
            throw std::domain_error("GlobalV::XC_Functional::get_func_type() =" + std::to_string(XC_Functional::get_func_type())
                + " unfinished in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
        }

        ModuleBase::timer::tick("PotHxcLR", "cal_v_eff");
    }

} // namespace elecstate
