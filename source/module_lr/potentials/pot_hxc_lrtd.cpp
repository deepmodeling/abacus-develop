#include "pot_hxc_lrtd.h"
#include "module_parameter/parameter.h"
#include "module_elecstate/potentials/H_Hartree_pw.h"
#include "module_base/timer.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include <set>
#include "module_lr/utils/lr_util.h"
#include "module_lr/utils/lr_util_xc.hpp"
#include "module_hamilt_pw/hamilt_pwdft/global.h"    // tmp, for pgrid
#define FXC_PARA_TYPE const double* const rho, ModuleBase::matrix& v_eff, const std::vector<int>& ispin_op = { 0,0 }
namespace LR
{
    // constructor for exchange-correlation kernel
    PotHxcLR::PotHxcLR(const std::string& xc_kernel, const ModulePW::PW_Basis& rho_basis, const UnitCell& ucell,
        const Charge& chg_gs/*ground state*/, const Parallel_Grid& pgrid,
        const SpinType& st, const std::vector<std::string>& lr_init_xc_kernel)
        :xc_kernel_(xc_kernel), tpiba_(ucell.tpiba), spin_type_(st), rho_basis_(rho_basis), nrxx_(chg_gs.nrxx),
        nspin_(PARAM.inp.nspin == 1 || (PARAM.inp.nspin == 4 && !PARAM.globalv.domag && !PARAM.globalv.domag_z) ? 1 : 2),
        pot_hartree_(LR_Util::make_unique<elecstate::PotHartree>(&rho_basis)),
        xc_kernel_components_(rho_basis, ucell, chg_gs, pgrid, nspin_, xc_kernel, lr_init_xc_kernel), //call XC_Functional::set_func_type and libxc
        xc_type_(XCType(XC_Functional::get_func_type()))
    {
        if (std::set<std::string>({ "lda", "pwlda", "pbe", "hse" }).count(xc_kernel)) { this->set_integral_func(this->spin_type_, this->xc_type_); }
    }

    void PotHxcLR::cal_v_eff(double** rho, const UnitCell& ucell, ModuleBase::matrix& v_eff, const std::vector<int>& ispin_op)
    {
        ModuleBase::TITLE("PotHxcLR", "cal_v_eff");
        ModuleBase::timer::tick("PotHxcLR", "cal_v_eff");
        auto& fxc = this->xc_kernel_components_;

        // Hartree
        switch (this->spin_type_)
        {
        case SpinType::S1: case SpinType::S2_updown:
            v_eff += elecstate::H_Hartree_pw::v_hartree(ucell, const_cast<ModulePW::PW_Basis*>(&this->rho_basis_), 1, rho);
            break;
        case SpinType::S2_singlet:
            v_eff += 2 * elecstate::H_Hartree_pw::v_hartree(ucell, const_cast<ModulePW::PW_Basis*>(&this->rho_basis_), 1, rho);
            break;
        default:
            break;
        }
        // XC
        if (this->xc_kernel_ == "rpa" || this->xc_kernel_ == "hf") { return; }    // no xc
#ifdef USE_LIBXC
        this->kernel_to_potential_[spin_type_](rho[0], v_eff, ispin_op);
#else
        throw std::domain_error("GlobalV::XC_Functional::get_func_type() =" + std::to_string(XC_Functional::get_func_type())
            + " unfinished in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
#endif
        ModuleBase::timer::tick("PotHxcLR", "cal_v_eff");
    }

    void PotHxcLR::set_integral_func(const SpinType& s, const XCType& xc)
    {
        auto& funcs = this->kernel_to_potential_;
        auto& fxc = this->xc_kernel_components_;
        if (xc == XCType::LDA) { switch (s)
        {
        case SpinType::S1:
            funcs[s] = [this, &fxc](FXC_PARA_TYPE)->void
                {
                    for (int ir = 0;ir < nrxx;++ir) { v_eff(0, ir) += ModuleBase::e2 * fxc.v2rho2.at(ir) * rho[ir]; }
                };
            break;
        case SpinType::S2_singlet:
            funcs[s] = [this, &fxc](FXC_PARA_TYPE)->void
                {
                    for (int ir = 0;ir < nrxx;++ir)
                    {
                        const int irs0 = 3 * ir;
                        const int irs1 = irs0 + 1;
                        v_eff(0, ir) += ModuleBase::e2 * (fxc.v2rho2.at(irs0) + fxc.v2rho2.at(irs1)) * rho[ir];
                    }
                };
            break;
        case SpinType::S2_triplet:
            funcs[s] = [this, &fxc](FXC_PARA_TYPE)->void
                {
                    for (int ir = 0;ir < nrxx;++ir)
                    {
                        const int irs0 = 3 * ir;
                        const int irs1 = irs0 + 1;
                        v_eff(0, ir) += ModuleBase::e2 * (fxc.v2rho2.at(irs0) - fxc.v2rho2.at(irs1)) * rho[ir];
                    }
                };
            break;
        case SpinType::S2_updown:
            funcs[s] = [this, &fxc](FXC_PARA_TYPE)->void
                {
                    assert(ispin_op.size() >= 2);
                    const int is = ispin_op[0] + ispin_op[1];
                    for (int ir = 0;ir < nrxx;++ir) { v_eff(0, ir) += ModuleBase::e2 * fxc.v2rho2.at(3 * ir + is) * rho[ir]; }
                };
            break;
        default:
            throw std::domain_error("SpinType =" + std::to_string(static_cast<int>(s))
                + " unfinished in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
            break;
        }
        } else if (xc == XCType::GGA || xc == XCType::HYB_GGA) { switch (s)
        {
        case SpinType::S1:
            funcs[s] = [this, &fxc](FXC_PARA_TYPE)->void
                {
                    // test: output drho
                    // double thr = 1e-1;
                    // auto out_thr = [this, &thr](const double* v) {
                    //     for (int ir = 0;ir < nrxx;++ir) if (std::abs(v[ir]) > thr) std::cout << v[ir] << " ";
                    //     std::cout << std::endl;};
                    // auto out_thr3 = [this, &thr](const std::vector<ModuleBase::Vector3<double>>& v) {
                    //     for (int ir = 0;ir < nrxx;++ir) if (std::abs(v.at(ir).x) > thr) std::cout << v.at(ir).x << " ";
                    //     std::cout << std::endl;
                    //     for (int ir = 0;ir < nrxx;++ir) if (std::abs(v.at(ir).y) > thr) std::cout << v.at(ir).y << " ";
                    //     std::cout << std::endl;
                    //     for (int ir = 0;ir < nrxx;++ir) if (std::abs(v.at(ir).z) > thr) std::cout << v.at(ir).z << " ";
                    //     std::cout << std::endl;};

                    std::vector<ModuleBase::Vector3<double>> drho(nrxx);    // transition density gradient
                    LR_Util::grad(rho, drho.data(), this->rho_basis_, this->tpiba_);

                    std::vector<double> vxc_tmp(nrxx, 0.0);

                    //1. $\partial E/\partial\rho = 2f^{\rho\sigma}*\nabla\rho*\rho_1+4f^{\sigma\sigma}\nabla\rho(\nabla\rho\cdot\nabla\rho_1)+2v^\sigma\nabla\rho_1$
                    std::vector<ModuleBase::Vector3<double>> e_drho(nrxx);
                    for (int ir = 0;ir < nrxx;++ir)
                    {
                        e_drho[ir] = -(fxc.v2rhosigma_2drho.at(ir) * rho[ir]
                            + fxc.v2sigma2_4drho.at(ir) * (fxc.drho_gs.at(ir) * drho.at(ir))
                            + drho.at(ir) * fxc.vsigma.at(ir) * 2.);
                    }
                    XC_Functional::grad_dot(e_drho.data(), vxc_tmp.data(), &this->rho_basis_, this->tpiba_);

                    // 2. $f^{\rho\rho}\rho_1+2f^{\rho\sigma}\nabla\rho\cdot\nabla\rho_1$
                    for (int ir = 0;ir < nrxx;++ir)
                    {
                        vxc_tmp[ir] += (fxc.v2rho2.at(ir) * rho[ir]
                            + fxc.v2rhosigma_2drho.at(ir) * drho.at(ir));
                    }
                    BlasConnector::axpy(nrxx, ModuleBase::e2, vxc_tmp.data(), 1, v_eff.c, 1);
                };
            break;
            // case SpinType::S2_singlet:
            //     break;
            // case SpinType::S2_triplet:
            //     break;
        default:
            throw std::domain_error("SpinType =" + std::to_string(static_cast<int>(s)) + "for GGA or HYB_GGA is unfinished in "
                + std::string(__FILE__) + " line " + std::to_string(__LINE__));
            break;
        }
        } else
        {
            throw std::domain_error("GlobalV::XC_Functional::get_func_type() =" + std::to_string(XC_Functional::get_func_type())
                + " unfinished in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
        }
    }
} // namespace LR
