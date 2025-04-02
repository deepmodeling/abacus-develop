#include "xc_kernel.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include "module_parameter/parameter.h"
#include "module_base/timer.h"
#include "module_lr/utils/lr_util.h"
#include "module_lr/utils/lr_util_xc.hpp"
#include <set>
#include <chrono>
#include "module_io/cube_io.h"
#ifdef USE_LIBXC
#include <xc.h>
#include "module_hamilt_general/module_xc/xc_functional_libxc.h"
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

LR::KernelXC::KernelXC(const ModulePW::PW_Basis& rho_basis,
    const UnitCell& ucell,
    const Charge& chg_gs,
    const Parallel_Grid& pgrid,
    const int& nspin,
    const std::string& kernel_name,
    const std::vector<std::string>& lr_init_xc_kernel,
    const bool openshell) :rho_basis_(rho_basis), openshell_(openshell)
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

template <typename T>
inline void add_op(const T* const src1, const T* const src2, T* const dst, const int size)
{
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096)
#endif
    for (int i = 0;i < size;++i)
    {
        dst[i] = src1[i] + src2[i];
    }
}
template <typename T>
inline void add_op(const std::vector<T>& src1, const std::vector<T>& src2, std::vector<T>& dst)
{
    assert(dst.size() >= src1.size() && src2.size() >= src1.size());
    add_op(src1.data(), src2.data(), dst.data(), src1.size());
}
template <typename T>
inline void add_assign_op(const std::vector<T>& src, std::vector<T>& dst)
{
    add_op(src, dst, dst);
}
template<typename Telement, typename Tscalar>
inline void cutoff_grid_data_spin2(std::vector<Telement>& func, const std::vector<Tscalar>& mask)
{
    const int& nrxx = mask.size() / 2;
    assert(func.size() % nrxx == 0 && func.size() / nrxx > 1);
    const int n_component = func.size() / nrxx;
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096)
#endif
    for (int ir = 0;ir < nrxx;++ir)
    {
        const int& i2 = 2 * ir;
        const int& istart = n_component * ir;
        std::for_each(func.begin() + istart, func.begin() + istart + n_component - 1, [&](Telement& f) { f *= mask[i2]; });    //spin-up
        std::for_each(func.begin() + istart + 1, func.begin() + istart + n_component, [&](Telement& f) { f *= mask[i2 + 1]; });    //spin-down
    }
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
            for (int ir = 0; ir < nrxx; ++ir) { rhor[ir] = rho[ir * nspin + is];
}
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
            for (int ir = 0; ir < nrxx; ++ir) {
                sigma[ir] = gradrho[0][ir] * gradrho[0][ir];
}
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
    this->vrho_.resize(nspin * nrxx, 0.);
    this->v2rho2_.resize(((1 == nspin) ? 1 : 3) * nrxx, 0.);//(nrxx* ((1 == nspin) ? 1 : 3)): 00, 01, 11
    if (is_gga)
    {
        this->vsigma_.resize(((1 == nspin) ? 1 : 3) * nrxx, 0.);//(nrxx*): 2 for rho * 3 for sigma: 00, 01, 02, 10, 11, 12
        this->v2rhosigma_.resize(((1 == nspin) ? 1 : 6) * nrxx, 0.); //(nrxx*): 2 for rho * 3 for sigma: 00, 01, 02, 10, 11, 12
        this->v2sigma2_.resize(((1 == nspin) ? 1 : 6) * nrxx, 0.);   //(nrxx* ((1 == nspin) ? 1 : 6)): 00, 01, 02, 11, 12, 22
    }
    //MetaGGA ...

    for (xc_func_type& func : funcs)
    {
        const double rho_threshold = 1E-6;
        const double grho_threshold = 1E-10;

        xc_func_set_dens_threshold(&func, rho_threshold);

        //cut off function
        const std::vector<double> sgn = XC_Functional_Libxc::cal_sgn(rho_threshold, grho_threshold, func, nspin, nrxx, rho, sigma);

        // Libxc interfaces overwrite (instead of add onto) the output arrays, so we need temporary copies
        std::vector<double> vrho_tmp(this->vrho_.size());
        std::vector<double> v2rho2_tmp(this->v2rho2_.size());
        std::vector<double> vsigma_tmp(this->vsigma_.size());
        std::vector<double> v2rhosigma_tmp(this->v2rhosigma_.size());
        std::vector<double> v2sigma2_tmp(this->v2sigma2_.size());
        switch (func.info->family)
        {
        case XC_FAMILY_LDA:
            xc_lda_vxc(&func, nrxx, rho.data(), vrho_tmp.data());
            xc_lda_fxc(&func, nrxx, rho.data(), v2rho2_tmp.data());
            break;
        case XC_FAMILY_GGA:
        case XC_FAMILY_HYB_GGA:
        {
            xc_gga_vxc(&func, nrxx, rho.data(), sigma.data(), vrho_tmp.data(), vsigma_tmp.data());
            xc_gga_fxc(&func, nrxx, rho.data(), sigma.data(), v2rho2_tmp.data(), v2rhosigma_tmp.data(), v2sigma2_tmp.data());
            // std::cout << "max element of v2sigma2_tmp: " << *std::max_element(v2sigma2_tmp.begin(), v2sigma2_tmp.end()) << std::endl;
            // std::cout << "rho corresponding to max element of v2sigma2_tmp: " << rho[(std::max_element(v2sigma2_tmp.begin(), v2sigma2_tmp.end()) - v2sigma2_tmp.begin()) / 6] << std::endl;
            // cut off by sgn
            cutoff_grid_data_spin2(vrho_tmp, sgn);
            cutoff_grid_data_spin2(vsigma_tmp, sgn);
            cutoff_grid_data_spin2(v2rho2_tmp, sgn);
            cutoff_grid_data_spin2(v2rhosigma_tmp, sgn);
            cutoff_grid_data_spin2(v2sigma2_tmp, sgn);
            break;
        }
        default:
            throw std::domain_error("func.info->family =" + std::to_string(func.info->family)
                + " unfinished in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
            break;
        }
        // add onto the total components
        // auto start = std::chrono::high_resolution_clock::now();
        add_assign_op(vrho_tmp, this->vrho_);
        add_assign_op(v2rho2_tmp, this->v2rho2_);
        add_assign_op(vsigma_tmp, this->vsigma_);
        add_assign_op(v2rhosigma_tmp, this->v2rhosigma_);
        add_assign_op(v2sigma2_tmp, this->v2sigma2_);
        // auto end = std::chrono::high_resolution_clock::now();
        // auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        // std::cout << "Time elapsed adding XC components: " << duration.count() << " ms\n";
    } // end for( xc_func_type &func : funcs )

    XC_Functional_Libxc::finish_func(funcs);

    // some postprocess if there're GGA funcs in the list
    if (is_gga)
    {
        const std::vector<double>& v2r2 = this->v2rho2_;
        const std::vector<double>& v2rs = this->v2rhosigma_;
        const std::vector<double>& v2s2 = this->v2sigma2_;
        const std::vector<double>& vs = this->vsigma_;
        const double tpiba2 = tpiba * tpiba;

        if (nspin == 1)
        {
            // 1. $2f^{\rho\sigma}*\nabla\rho$
            this->v2rhosigma_2drho_.resize(nrxx);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096)
#endif
            for (size_t i = 0; i < nrxx; ++i)
            {
                this->v2rhosigma_2drho_[i] = gradrho[0][i] * v2rs[i] * 2.;
            }

            // 2. $4f^{\sigma\sigma}*\nabla\rho$
            this->v2sigma2_4drho_.resize(nrxx);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096)
#endif
            for (size_t i = 0; i < nrxx; ++i)
            {
                this->v2sigma2_4drho_[i] = gradrho[0][i] * v2s2[i] * 4.;
            }
        }
        else if (2 == nspin)    //close-shell
        {
            if (!openshell_)
            {
                this->v2rhosigma_drho_singlet_.resize(nrxx);
                this->v2sigma2_drho_singlet_.resize(nrxx);
                this->v2rhosigma_drho_triplet_.resize(nrxx);
                this->v2sigma2_drho_triplet_.resize(nrxx);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096)
#endif
                for (int i = 0;i < nrxx;++i)
                {
                    const int istart = i * 6;
                    this->v2rhosigma_drho_singlet_[i] = gradrho[0][i] * (v2rs[istart] + v2rs[istart + 1] + v2rs[istart + 2]) * 2.;
                    this->v2sigma2_drho_singlet_[i] = gradrho[0][i] * (v2s2[istart] * 2. + v2s2[istart + 1] * 3. + v2s2[istart + 2] * 2. + v2s2[istart + 3] + v2s2[istart + 4]) * 2.;
                    this->v2rhosigma_drho_triplet_[i] = gradrho[0][i] * (v2rs[istart] - v2rs[istart + 2]) * 2.;
                    this->v2sigma2_drho_triplet_[i] = gradrho[0][i] * (v2s2[istart] * 2. + v2s2[istart + 1] - v2s2[istart + 2] * 2. - v2s2[istart + 4]) * 2.;
                }
            }
            else
            {
                this->v2rhosigma_drho_uu_.resize(nrxx);
                this->v2rhosigma_drho_ud_.resize(nrxx);
                this->v2rhosigma_drho_du_.resize(nrxx);
                this->v2rhosigma_drho_dd_.resize(nrxx);
                this->v2sigma2_drho_uu_u_.resize(nrxx);
                this->v2sigma2_drho_uu_d_.resize(nrxx);
                this->v2sigma2_drho_ud_u_.resize(nrxx);
                this->v2sigma2_drho_ud_d_.resize(nrxx);
                this->v2sigma2_drho_du_u_.resize(nrxx);
                this->v2sigma2_drho_du_d_.resize(nrxx);
                this->v2sigma2_drho_dd_u_.resize(nrxx);
                this->v2sigma2_drho_dd_d_.resize(nrxx);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096)
#endif
                for (int i = 0;i < nrxx;++i)
                {
                    const int istart = i * 6;
                    this->v2rhosigma_drho_uu_[i] = gradrho[0][i] * v2rs[istart] * 2. + gradrho[1][i] * v2rs[istart + 1];
                    this->v2rhosigma_drho_ud_[i] = gradrho[0][i] * v2rs[istart + 1] + gradrho[1][i] * v2rs[istart + 2] * 2.;
                    this->v2rhosigma_drho_du_[i] = gradrho[0][i] * v2rs[istart + 3] * 2. + gradrho[1][i] * v2rs[istart + 4];
                    this->v2rhosigma_drho_dd_[i] = gradrho[0][i] * v2rs[istart + 4] + gradrho[1][i] * v2rs[istart + 5] * 2.;
                    this->v2sigma2_drho_uu_u_[i] = gradrho[0][i] * v2s2[istart] * 4. + gradrho[1][i] * v2s2[istart + 1] * 2.;
                    this->v2sigma2_drho_uu_d_[i] = gradrho[0][i] * v2s2[istart + 1] * 2. + gradrho[1][i] * v2s2[istart + 3];
                    this->v2sigma2_drho_ud_u_[i] = gradrho[0][i] * v2s2[istart + 1] * 2. + gradrho[1][i] * v2s2[istart + 3];
                    this->v2sigma2_drho_ud_d_[i] = gradrho[0][i] * v2s2[istart + 2] * 4. + gradrho[1][i] * v2s2[istart + 4] * 2.;
                    this->v2sigma2_drho_du_u_[i] = gradrho[1][i] * v2s2[istart + 2] * 4. + gradrho[0][i] * v2s2[istart + 1] * 2.;
                    this->v2sigma2_drho_du_d_[i] = gradrho[1][i] * v2s2[istart + 4] * 2. + gradrho[0][i] * v2s2[istart + 3];
                    this->v2sigma2_drho_dd_u_[i] = gradrho[1][i] * v2s2[istart + 4] * 2. + gradrho[0][i] * v2s2[istart + 3];
                    this->v2sigma2_drho_dd_d_[i] = gradrho[1][i] * v2s2[istart + 5] * 4. + gradrho[0][i] * v2s2[istart + 4] * 2.;
                }
            }
        }
        else
        {
            throw std::domain_error("nspin =" + std::to_string(nspin)
                + " unfinished in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
        }
        this->drho_gs_ = std::move(gradrho);
    }
    if (PARAM.inp.nspin == 1 || PARAM.inp.nspin == 2) {
        return;
    // else if (4 == PARAM.inp.nspin)
    } else//NSPIN != 1,2,4 is not supported
    {
        throw std::domain_error("PARAM.inp.nspin =" + std::to_string(PARAM.inp.nspin)
            + " unfinished in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
    }
}
#endif