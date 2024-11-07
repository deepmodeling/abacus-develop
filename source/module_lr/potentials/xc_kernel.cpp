#include "xc_kernel.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include "module_parameter/parameter.h"
#include "module_base/timer.h"
#include "module_lr/utils/lr_util.h"
#include "module_lr/utils/lr_util_xc.hpp"
#include <set>
#include "module_io/cube_io.h"
#ifdef USE_LIBXC
#include <xc.h>
#include "module_hamilt_general/module_xc/xc_functional_libxc.h"
#endif

LR::KernelXC::KernelXC(const ModulePW::PW_Basis& rho_basis,
    const UnitCell& ucell,
    const Charge& chg_gs,
    const Parallel_Grid& pgrid,
    const int& nspin,
    const std::string& kernel_name,
    const std::vector<std::string>& lr_init_xc_kernel) :rho_basis_(rho_basis)
{
    if (!std::set<std::string>({ "lda", "pwlda", "pbe", "hse" }).count(kernel_name)) { return; }
    XC_Functional::set_xc_type(kernel_name);    // for hse, (1-alpha) and omega are set here

    const int& nrxx = rho_basis.nrxx;
    if (lr_init_xc_kernel[0] == "file")
    {
        const std::set<std::string> lda_xc = { "lda", "pwlda" };
        assert(lda_xc.count(kernel_name));
        const int n_component = (1 == nspin) ? 1 : 3;   // spin components of fxc: (uu, ud=du, dd) when nspin=2
        this->v2rho2_.resize(n_component * nrxx);
        // read fxc adn add to xc_kernel_components
        assert(lr_init_xc_kernel.size() >= n_component + 1);
        for (int is = 0;is < n_component;++is)
        {
            double ef = 0.0;
            int prenspin = 1;
            std::vector<double> v2rho2_tmp(nrxx);
            ModuleIO::read_vdata_palgrid(pgrid, GlobalV::MY_RANK, GlobalV::ofs_running, lr_init_xc_kernel[is + 1],
                v2rho2_tmp.data(), ucell.nat);
            for (int ir = 0;ir < nrxx;++ir) { this->v2rho2_[ir * n_component + is] = v2rho2_tmp[ir]; }
        }
        return;
    }

#ifdef USE_LIBXC
    if (lr_init_xc_kernel[0] == "from_charge_file")
    {
        assert(lr_init_xc_kernel.size() >= 2);
        double** rho_for_fxc;
        LR_Util::_allocate_2order_nested_ptr(rho_for_fxc, nspin, nrxx);
        double ef = 0.0;
        int prenspin = 1;
        for (int is = 0;is < nspin;++is)
        {
            const std::string file = lr_init_xc_kernel[lr_init_xc_kernel.size() > nspin ? 1 + is : 1];
            ModuleIO::read_vdata_palgrid(pgrid, GlobalV::MY_RANK, GlobalV::ofs_running, file,
                rho_for_fxc[is], ucell.nat);
        }
        this->f_xc_libxc(nspin, ucell.omega, ucell.tpiba, rho_for_fxc, chg_gs.rho_core);
        LR_Util::_deallocate_2order_nested_ptr(rho_for_fxc, nspin);
    }
    else { this->f_xc_libxc(nspin, ucell.omega, ucell.tpiba, chg_gs.rho, chg_gs.rho_core); }
#else 
    ModuleBase::WARNING_QUIT("KernelXC", "to calculate xc-kernel in LR-TDDFT, compile with LIBXC");
#endif
}

#ifdef USE_LIBXC
void LR::KernelXC::f_xc_libxc(const int& nspin, const double& omega, const double& tpiba, const double* const* const rho_gs, const double* const rho_core)
{
    ModuleBase::TITLE("XC_Functional", "f_xc_libxc");
    ModuleBase::timer::tick("XC_Functional", "f_xc_libxc");
    // https://www.tddft.org/programs/libxc/manual/libxc-5.1.x/

    assert(nspin == 1 || nspin == 2);

    std::vector<xc_func_type> funcs = XC_Functional_Libxc::init_func(
        XC_Functional::get_func_id(),
        (1 == nspin) ? XC_UNPOLARIZED : XC_POLARIZED);
    const int& nrxx = rho_basis_.nrxx;

    // converting rho (extract it as a subfuntion in the future)
    // -----------------------------------------------------------------------------------
    std::vector<double> rho(nspin * nrxx);    // r major / spin contigous

#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 1024)
#endif
    for (int is = 0; is < nspin; ++is) { for (int ir = 0; ir < nrxx; ++ir) { rho[ir * nspin + is] = rho_gs[is][ir]; } }
    if (rho_core)
    {
        const double fac = 1.0 / nspin;
        for (int is = 0; is < nspin; ++is) { for (int ir = 0; ir < nrxx; ++ir) { rho[ir * nspin + is] += fac * rho_core[ir]; } }
    }

    // -----------------------------------------------------------------------------------
    // for GGA
    const bool is_gga = std::any_of(funcs.begin(), funcs.end(), [](const xc_func_type& f) { return f.info->family == XC_FAMILY_GGA || f.info->family == XC_FAMILY_HYB_GGA; });

    std::vector<std::vector<ModuleBase::Vector3<double>>> gradrho;  // \nabla \rho
    std::vector<double> sigma;  // |\nabla\rho|^2
    std::vector<double> sgn;        // sgn for threshold mask
    if (is_gga)
    {
        // 0. set up sgn for threshold mask
        // in the case of GGA correlation for polarized case,
        // a cutoff for grho is required to ensure that libxc gives reasonable results

        // 1. \nabla \rho
        gradrho.resize(nspin);
        for (int is = 0; is < nspin; ++is)
        {
            std::vector<double> rhor(nrxx);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
            for (int ir = 0; ir < nrxx; ++ir) rhor[ir] = rho[ir * nspin + is];
            gradrho[is].resize(nrxx);
            LR_Util::grad(rhor.data(), gradrho[is].data(), rho_basis_, tpiba);
        }
        // 2. |\nabla\rho|^2
        sigma.resize(nrxx * ((1 == nspin) ? 1 : 3));
        if (1 == nspin)
        {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
            for (int ir = 0; ir < nrxx; ++ir)
                sigma[ir] = gradrho[0][ir] * gradrho[0][ir];
        }
        else
        {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 256)
#endif
            for (int ir = 0; ir < nrxx; ++ir)
            {
                sigma[ir * 3] = gradrho[0][ir] * gradrho[0][ir];
                sigma[ir * 3 + 1] = gradrho[0][ir] * gradrho[1][ir];
                sigma[ir * 3 + 2] = gradrho[1][ir] * gradrho[1][ir];
            }
        }
    }
    // -----------------------------------------------------------------------------------
    //==================== XC Kernels (f_xc)=============================
    this->vrho_.resize(nspin * nrxx);
    this->v2rho2_.resize(((1 == nspin) ? 1 : 3) * nrxx);//(nrxx* ((1 == nspin) ? 1 : 3)): 00, 01, 11
    if (is_gga)
    {
        this->vsigma_.resize(((1 == nspin) ? 1 : 3) * nrxx);//(nrxx*): 2 for rho * 3 for sigma: 00, 01, 02, 10, 11, 12
        this->v2rhosigma_.resize(((1 == nspin) ? 1 : 6) * nrxx); //(nrxx*): 2 for rho * 3 for sigma: 00, 01, 02, 10, 11, 12
        this->v2sigma2_.resize(((1 == nspin) ? 1 : 6) * nrxx);   //(nrxx* ((1 == nspin) ? 1 : 6)): 00, 01, 02, 11, 12, 22
    }
    //MetaGGA ...

    for (xc_func_type& func : funcs)
    {
        const double rho_threshold = 1E-6;
        const double grho_threshold = 1E-10;

        xc_func_set_dens_threshold(&func, rho_threshold);

        //cut off grho if not LDA (future subfunc)

        switch (func.info->family)
        {
        case XC_FAMILY_LDA:
            xc_lda_vxc(&func, nrxx, rho.data(), vrho_.data());
            xc_lda_fxc(&func, nrxx, rho.data(), v2rho2_.data());
            break;
        case XC_FAMILY_GGA:
        case XC_FAMILY_HYB_GGA:
            xc_gga_vxc(&func, nrxx, rho.data(), sigma.data(), vrho_.data(), vsigma_.data());
            xc_gga_fxc(&func, nrxx, rho.data(), sigma.data(), v2rho2_.data(), v2rhosigma_.data(), v2sigma2_.data());
            break;
        default:
            throw std::domain_error("func.info->family =" + std::to_string(func.info->family)
                + " unfinished in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
            break;
        }
        // some formulas for GGA
        if (func.info->family == XC_FAMILY_GGA || func.info->family == XC_FAMILY_HYB_GGA)
        {
            const std::vector<double>& v2r2 = this->v2rho2_;
            const std::vector<double>& v2rs = this->v2rhosigma_;
            const std::vector<double>& v2s2 = this->v2sigma2_;
            const std::vector<double>& vs = this->vsigma_;
            const double tpiba2 = tpiba * tpiba;

            if (1 == nspin)
            {
                using V3 = ModuleBase::Vector3<double>;
                // 0. drho
                this->drho_gs_ = gradrho[0];
                // 1. $2f^{\rho\sigma}*\nabla\rho$
                this->v2rhosigma_2drho_.resize(nrxx);
                std::transform(gradrho[0].begin(), gradrho[0].end(), v2rs.begin(), this->v2rhosigma_2drho_.begin(),
                    [](const V3& a, const V3& b) {return a * b * 2.; });
                // 2. $4f^{\sigma\sigma}*\nabla\rho$
                this->v2sigma2_4drho_.resize(nrxx);
                std::transform(sigma.begin(), sigma.end(), v2s2.begin(), v2sigma2_4drho_.begin(),
                    [](const V3& a, const V3& b) {return a * b * 4.; });
            }
            // else if (2 == nspin)    // wrong, to be fixed
            // {
            //     // 1. $\nabla\cdot(f^{\rho\sigma}*\nabla\rho)$
            //     std::vector<double> div_v2rhosigma_gdrho_r(3 * nrxx);
            //     std::vector<ModuleBase::Vector3<double>> v2rhosigma_gdrho_r(3 * nrxx);
            //     for (int ir = 0; ir < nrxx; ++ir)
            //     {
            //         v2rhosigma_gdrho_r[ir] = gradrho[0][ir] * v2rs.at(ir * 6) * 4.0
            //             + gradrho[1][ir] * v2rs.at(ir * 6 + 1) * 2.0;   //up-up
            //         v2rhosigma_gdrho_r[nrxx + ir] = gradrho[0][ir] * (v2rs.at(ir * 6 + 3) * 2.0 + v2rs.at(ir * 6 + 1))
            //             + gradrho[1][ir] * (v2rs.at(ir * 6 + 2) * 2.0 + v2rs.at(ir * 6 + 4));   //up-down
            //         v2rhosigma_gdrho_r[2 * nrxx + ir] = gradrho[1][ir] * v2rs.at(ir * 6 + 5) * 4.0
            //             + gradrho[0][ir] * v2rs.at(ir * 6 + 4) * 2.0;   //down-down
            //     }
            //     for (int isig = 0;isig < 3;++isig)
            //         XC_Functional::grad_dot(v2rhosigma_gdrho_r.data() + isig * nrxx, div_v2rhosigma_gdrho_r.data() + isig * nrxx, chg_gs.rhopw, tpiba);
            //     // 2. $\nabla^2(f^{\sigma\sigma}*\sigma)$
            //     std::vector<double> v2sigma2_sigma_r(3 * nrxx);
            //     for (int ir = 0; ir < nrxx; ++ir)
            //     {
            //         v2sigma2_sigma_r[ir] = v2s2.at(ir * 6) * sigma[ir * 3] * 4.0
            //             + v2s2.at(ir * 6 + 1) * sigma[ir * 3 + 1] * 4.0
            //             + v2s2.at(ir * 6 + 3) * sigma[ir * 3 + 2];   //up-up
            //         v2sigma2_sigma_r[nrxx + ir] = v2s2.at(ir * 6 + 1) * sigma[ir * 3] * 2.0
            //             + v2s2.at(ir * 6 + 4) * sigma[ir * 3 + 2] * 2.0
            //             + (v2s2.at(ir * 6 + 2) * 4.0 + v2s2.at(ir * 6 + 3)) * sigma[ir * 3 + 1];   //up-down
            //         v2sigma2_sigma_r[2 * nrxx + ir] = v2s2.at(ir * 6 + 5) * sigma[ir * 3 + 2] * 4.0
            //             + v2s2.at(ir * 6 + 4) * sigma[ir * 3 + 1] * 4.0
            //             + v2s2.at(ir * 6 + 3) * sigma[ir * 3];   //down-down
            //     }
            //     for (int isig = 0;isig < 3;++isig)
            //         LR_Util::lapl(v2sigma2_sigma_r.data() + isig * nrxx, v2sigma2_sigma_r.data() + isig * nrxx, *(chg_gs.rhopw), tpiba2);
            //     // 3. $\nabla^2(v^\sigma)$
            //     std::vector<double> lap_vsigma(3 * nrxx);
            //     for (int ir = 0;ir < nrxx;++ir)
            //     {
            //         lap_vsigma[ir] = vs.at(ir * 3) * 2.0;
            //         lap_vsigma[nrxx + ir] = vs.at(ir * 3 + 1);
            //         lap_vsigma[2 * nrxx + ir] = vs.at(ir * 3 + 2) * 2.0;
            //     }
            //     for (int isig = 0;isig < 3;++isig)
            //         LR_Util::lapl(lap_vsigma.data() + isig * nrxx, lap_vsigma.data() + isig * nrxx, *(chg_gs.rhopw), tpiba2);
            //     // add to v2rho2_
            //     BlasConnector::axpy(3 * nrxx, 1.0, v2r2.data(), 1, to_mul_rho_.data(), 1);
            //     BlasConnector::axpy(3 * nrxx, -1.0, div_v2rhosigma_gdrho_r.data(), 1, to_mul_rho_.data(), 1);
            //     BlasConnector::axpy(3 * nrxx, 1.0, v2sigma2_sigma_r.data(), 1, to_mul_rho_.data(), 1);
            //     BlasConnector::axpy(3 * nrxx, 1.0, lap_vsigma.data(), 1, to_mul_rho_.data(), 1);
            // }
            else
            {
                throw std::domain_error("nspin =" + std::to_string(nspin)
                    + " unfinished in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
            }
        }
    } // end for( xc_func_type &func : funcs )
    XC_Functional_Libxc::finish_func(funcs);

    if (1 == PARAM.inp.nspin || 2 == PARAM.inp.nspin) return;
    // else if (4 == PARAM.inp.nspin)
    else//NSPIN != 1,2,4 is not supported
    {
        throw std::domain_error("PARAM.inp.nspin =" + std::to_string(PARAM.inp.nspin)
            + " unfinished in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
    }
}
#endif