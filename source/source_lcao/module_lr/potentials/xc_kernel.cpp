#include "xc_kernel.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_io/module_parameter/parameter.h"
#include "source_base/timer.h"
#include "source_lcao/module_lr/utils/lr_util.h"
#include "source_lcao/module_lr/utils/lr_util_xc.hpp"
#include <set>
#include <chrono>
#include "source_io/module_output/cube_io.h"
#ifdef __LIBXC
#include <xc.h>
#include "source_hamilt/module_xc/libxc_abacus.h"
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
    const bool openshell,
    const int gxc_spin) :rho_basis_(rho_basis), openshell_(openshell), gxc_spin_(gxc_spin)
{
    if (!LR_Util::has_local_xc(kernel_name)) { return; }
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

#ifdef __LIBXC
    if (lr_init_xc_kernel[0] == "from_charge_file")
    {
        assert(lr_init_xc_kernel.size() >= 2);
        double** rho_for_fxc = nullptr;
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

#ifdef __LIBXC
void LR::KernelXC::f_xc_libxc(const int& nspin, const double& omega, const double& tpiba, const double* const* const rho_gs, const double* const rho_core)
{
    ModuleBase::TITLE("XC_Functional", "f_xc_libxc");
    ModuleBase::timer::start("XC_Functional", "f_xc_libxc");
    // https://www.tddft.org/programs/libxc/manual/libxc-5.1.x/

    assert(nspin == 1 || nspin == 2);

    double hybrid_alpha = 0.0;
    double hse_omega = 0.0;
    std::vector<xc_func_type> funcs = XC_Functional_Libxc::init_func(
        XC_Functional::get_func_id(),
        (1 == nspin) ? XC_UNPOLARIZED : XC_POLARIZED,
        hybrid_alpha,
        hse_omega);
    const int& nrxx = rho_basis_.nrxx;
    const bool is_gga = std::any_of(funcs.begin(), funcs.end(), [](const xc_func_type& f) { return f.info->family == XC_FAMILY_GGA || f.info->family == XC_FAMILY_HYB_GGA; });
    // The third-order kernel exists for exactly one purpose: building the $g^{xc}$ coefficients
    // below. If none were requested, skip it. Openshell waits for future implementation.
    const bool need_kxc = (this->gxc_spin_ != GxcSpin::NoGxc) && !this->openshell_;

    std::vector<double> rho(nspin * nrxx);    // r major / spin contigous
    // for GGA
    std::vector<std::vector<ModuleBase::Vector3<double>>> gradrho;  // \nabla \rho
    std::vector<double> sigma;  // |\nabla\rho|^2
    std::vector<double> sgn;        // sgn for threshold mask
    this->get_rho_drho_sigma(nspin, tpiba, rho_gs, rho_core, is_gga, rho, gradrho, sigma);

    //==================== XC Kernels (f_xc)=============================
    this->vrho_.resize(nspin * nrxx, 0.);
    this->v2rho2_.resize(((1 == nspin) ? 1 : 3) * nrxx, 0.);//(nrxx* ((1 == nspin) ? 1 : 3)): 00, 01, 11
    if (need_kxc)
    {
        this->v3rho3_.resize(((1 == nspin) ? 1 : 4) * nrxx, 0.);//(nrxx* ((1 == nspin) ? 1 : 4)): 000, 001, 011, 111
    }
    if (is_gga)
    {
        this->vsigma_.resize(((1 == nspin) ? 1 : 3) * nrxx, 0.);//(nrxx*): 2 for rho * 3 for sigma: 00, 01, 02, 10, 11, 12
        this->v2rhosigma_.resize(((1 == nspin) ? 1 : 6) * nrxx, 0.); //(nrxx*): 2 for rho * 3 for sigma: 00, 01, 02, 10, 11, 12
        this->v2sigma2_.resize(((1 == nspin) ? 1 : 6) * nrxx, 0.);   //(nrxx* ((1 == nspin) ? 1 : 6)): 00, 01, 02, 11, 12, 22
        if (need_kxc)
        {
            this->v3rho2sigma_.resize(((1 == nspin) ? 1 : 9) * nrxx, 0.); //000, 001, 002, 010, 011, 012, 110, 111, 112
            this->v3rhosigma2_.resize(((1 == nspin) ? 1 : 12) * nrxx, 0.);   //000, 001, 002, 011, 012, 022, 100, 101, 102, 111, 112, 122
            this->v3sigma3_.resize(((1 == nspin) ? 1 : 10) * nrxx, 0.);//000, 001, 002, 011, 012, 022, 111, 112, 122, 222
        }
    }
    //MetaGGA ...

    for (xc_func_type& func : funcs)
    {
        // These used to be 1E-6 / 1E-10, the values the ground-state SCF uses. That is far too
        // aggressive for LR *gradients*: `cal_sgn` zeroes the kernel wherever rho < rho_threshold,
        // and in a molecule-in-a-big-box most of the grid is below 1E-6, while the exact derivative
        // the finite difference measures has no such truncation.
        //
        // Measured on H2/DZP TDRPA@LDA, where K^T = 0 makes the triplet gradient identical to
        // d(eps_a - eps_i)/dx and the reference is therefore exact to ~13 digits:
        //     thresholds 1E-6 /1E-10 : max |err| = 0.0159 eV/Ang over 9 states
        //     thresholds 1E-14/1E-20 : max |err| = 0.0002 eV/Ang  (at the force printout precision)
        // and on H2/DZP TDLDA over 18 states: 0.038 -> 0.00095 eV/Ang (the latter is the finite
        // difference's own resolution). H2/SZ was insensitive either way -- its compact 1s-only
        // basis puts no T+D^Z density in the truncated region, which is why the problem only shows
        // up once diffuse/p functions enter.
        //
        // CAVEAT: only LDA has been checked at these thresholds. The 1E-6/1E-10 pair exists because
        // GGA *correlation* can misbehave at very low density; if PBE turns out to need protection,
        // the fix is a separate threshold for the gradient path, not a return to 1E-6 everywhere.
        const double rho_threshold = 1E-14;
        const double grho_threshold = 1E-20;

        xc_func_set_dens_threshold(&func, rho_threshold);

        //cut off function
        const std::vector<double> sgn = XC_Functional_Libxc::cal_sgn(rho_threshold, grho_threshold, func, nspin, nrxx, rho, sigma);

        // Libxc interfaces overwrite (instead of add onto) the output arrays, so we need temporary copies
        std::vector<double> vrho_tmp(this->vrho_.size());
        std::vector<double> v2rho2_tmp(this->v2rho2_.size());
        std::vector<double> vsigma_tmp(this->vsigma_.size());
        std::vector<double> v2rhosigma_tmp(this->v2rhosigma_.size());
        std::vector<double> v2sigma2_tmp(this->v2sigma2_.size());
        // The kxc arrays used to be handed to libxc directly. Since the libxc interfaces *overwrite*
        // their output (that is why every other component already goes through a temporary), only the
        // last functional of a composite survived -- e.g. for PBE = XC_GGA_X_PBE + XC_GGA_C_PBE the
        // exchange part of the third derivative was silently dropped. They are accumulated now too.
        // Note these are left uncut by `sgn`: for nspin=1 the cutoff is a no-op anyway (a single
        // spin component), and for nspin=2 `cutoff_grid_data_spin2` does not match the component
        // layout of the third derivatives.
        std::vector<double> v3rho3_tmp(this->v3rho3_.size());
        std::vector<double> v3rho2sigma_tmp(this->v3rho2sigma_.size());
        std::vector<double> v3rhosigma2_tmp(this->v3rhosigma2_.size());
        std::vector<double> v3sigma3_tmp(this->v3sigma3_.size());
        switch (func.info->family)
        {
        case XC_FAMILY_LDA:
            xc_lda_vxc(&func, nrxx, rho.data(), vrho_tmp.data());
            xc_lda_fxc(&func, nrxx, rho.data(), v2rho2_tmp.data());
            if (need_kxc)
            {
                xc_lda_kxc(&func, nrxx, rho.data(), v3rho3_tmp.data());
            }
            break;
        case XC_FAMILY_GGA:
        case XC_FAMILY_HYB_GGA:
        {
            xc_gga_vxc(&func, nrxx, rho.data(), sigma.data(), vrho_tmp.data(), vsigma_tmp.data());
            xc_gga_fxc(&func, nrxx, rho.data(), sigma.data(), v2rho2_tmp.data(), v2rhosigma_tmp.data(), v2sigma2_tmp.data());
            // std::cout << "max element of v2sigma2_tmp: " << *std::max_element(v2sigma2_tmp.begin(), v2sigma2_tmp.end()) << std::endl;
            // std::cout << "rho corresponding to max element of v2sigma2_tmp: " << rho[(std::max_element(v2sigma2_tmp.begin(), v2sigma2_tmp.end()) - v2sigma2_tmp.begin()) / 6] << std::endl;
            // cut off by sgn. nspin=2 only: `cutoff_grid_data_spin2` assumes >1 component per
            // grid point (it asserts on it), and at nspin=1 there is exactly one, for which both
            // of its `for_each` ranges are empty -- the cutoff is a no-op anyway. 
            if (nspin == 2)
            {
                cutoff_grid_data_spin2(vrho_tmp, sgn);
                cutoff_grid_data_spin2(vsigma_tmp, sgn);
                cutoff_grid_data_spin2(v2rho2_tmp, sgn);
                cutoff_grid_data_spin2(v2rhosigma_tmp, sgn);
                cutoff_grid_data_spin2(v2sigma2_tmp, sgn);
            }
            if (need_kxc)
            {
                xc_gga_kxc(&func, nrxx, rho.data(), sigma.data(),
                    v3rho3_tmp.data(),
                    v3rho2sigma_tmp.data(),
                    v3rhosigma2_tmp.data(),
                    v3sigma3_tmp.data());
            }
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
        if (need_kxc)
        {
            add_assign_op(v3rho3_tmp, this->v3rho3_);
            add_assign_op(v3rho2sigma_tmp, this->v3rho2sigma_);
            add_assign_op(v3rhosigma2_tmp, this->v3rhosigma2_);
            add_assign_op(v3sigma3_tmp, this->v3sigma3_);
        }
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

            // The third-order kernels contracted with $\nabla\rho$ used to be pre-built here;
            // they are now folded into `GxcCoef::e_*` (with $\nabla\rho$ factored back out)
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
    }

    // Build the $v^{(2)}$ coefficient sets. Must happen before `gradrho` is moved away below.
    // Only the combinations the caller asked for: at nspin=2 each set costs 10 doubles per grid
    // point for a GGA, so building an unused one is a third of the whole third-order kernel.
    this->nspin_ = nspin;
    if (!openshell_ && this->gxc_spin_ != GxcSpin::NoGxc)
    {
        // No `gradrho` here: the divergence coefficients are stored with $\nabla\rho$ factored
        // out (see `GxcCoef`), so the spin algebra below is purely local.
        // nspin=1 has no triplet: `gxc()` returns `gxc_s_` whatever is asked, so build it for any request.
        if (nspin == 1 || (this->gxc_spin_ & GxcSpin::Singlet)) // Singlet or Bothspin
        {
            this->build_gxc_coef(this->gxc_s_, /*triplet=*/false, nspin, is_gga);
        }
        if (nspin == 2 && (this->gxc_spin_ & GxcSpin::Triplet)) // Triplet or Bothspin
        {
            this->build_gxc_coef(this->gxc_t_, /*triplet=*/true, nspin, is_gga);
        }
    }
    if (is_gga) { this->drho_gs_ = std::move(gradrho); }
    ModuleBase::timer::end("XC_Functional", "f_xc_libxc");
}

void LR::KernelXC::get_rho_drho_sigma(const int& nspin,
    const double& tpiba,
    const double* const* const rho_gs,
    const double* const rho_core,
    const bool& is_gga,
    std::vector<double>& rho,
    std::vector<std::vector<ModuleBase::Vector3<double>>>& gradrho,
    std::vector<double>& sigma)
{
    const int nrxx = rho_basis_.nrxx;
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 1024)
#endif
    for (int is = 0; is < nspin; ++is) { for (int ir = 0; ir < nrxx; ++ir) { rho[ir * nspin + is] = rho_gs[is][ir]; } }
    if (rho_core)
    {
        const double fac = 1.0 / nspin;
        for (int is = 0; is < nspin; ++is) { for (int ir = 0; ir < nrxx; ++ir) { rho[ir * nspin + is] += fac * rho_core[ir]; } }
    }
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
            for (int ir = 0; ir < nrxx; ++ir) {
                rhor[ir] = rho[ir * nspin + is];
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
}


// ---------------------------------------------------------------------------------------------
// The spin algebra of $v^{(2)}$, for both the singlet and the triplet combination.
//
// Both combinations come from the SAME lambda-expansion; they differ only in five weight vectors,
// which is why they share this one routine:
//
//   perturbation   singlet: rho_u += L*s, rho_d += L*s     triplet: rho_u += L*s, rho_d -= L*s
//   eta_sigma      (+1, +1)                                (+1, -1)     d(rho_sigma)/dL / s
//   mu_{ab}        (+1, +1, +1)                            (+1,  0, -1) d(sigma_ab)/dL / (2t)
//   nu_{ab}        (+1, +1, +1)                            (+1, -1, +1) d2(sigma_ab)/dL2 / (2q)
//   theta_{ab}     ( 2,  1,  0)                            ( 2,  1,  0) from 2*e_{s_uu}*grad rho_u
//   theta~_{ab}    ( 2,  1,  0)                            ( 2, -1,  0)   + e_{s_ud}*grad rho_d
//
// theta~ differs from theta only for the triplet, and only because grad rho_d there carries
// d(grad rho_d)/dL = -grad rho^1: W'' picks up 4u'_uu - 2u'_ud where W' picks up 2u'_uu + u'_ud.
// In the singlet the two coincide, which is why the singlet formula looked tidier than it is.
// Getting theta~ wrong is invisible in every singlet test.
namespace
{
    // libxc component indices for nspin=2.
    // sigma types: 0=uu, 1=ud, 2=dd. Unordered sigma pairs -> v2sigma2 / (per-rho block of) v3rhosigma2:
    constexpr int p2[3][3] = { {0,1,2},{1,3,4},{2,4,5} };
    // Unordered sigma triples -> v3sigma3: (000)(001)(002)(011)(012)(022)(111)(112)(122)(222)
    constexpr int p3[3][3][3] = {
        { {0,1,2},{1,3,4},{2,4,5} },
        { {1,3,4},{3,6,7},{4,7,8} },
        { {2,4,5},{4,7,8},{5,8,9} } };
}

void LR::KernelXC::build_gxc_coef(GxcCoef& dst, const bool triplet, const int& nspin, const bool& is_gga)
{
    const int& nrxx = rho_basis_.nrxx;
    const double eta[2] = { 1., triplet ? -1. : 1. };
    const double mu[3] = { 1., triplet ? 0. : 1., triplet ? -1. : 1. };
    const double nu[3] = { 1., triplet ? -1. : 1., 1. };
    const double th[3] = { 2., 1., 0. };
    const double tht[3] = { 2., triplet ? -1. : 1., 0. };

    dst.a_s2.resize(nrxx, 0.);
    if (is_gga)
    {
        dst.a_st.resize(nrxx, 0.); dst.a_t2.resize(nrxx, 0.); dst.a_q.resize(nrxx, 0.);
        dst.c_s.resize(nrxx, 0.); dst.c_t.resize(nrxx, 0.);
        dst.e_s2.resize(nrxx); dst.e_st.resize(nrxx); dst.e_t2.resize(nrxx); dst.e_q.resize(nrxx);
    }
    const std::vector<double>& v2rs = this->v2rhosigma_;
    const std::vector<double>& v2s2 = this->v2sigma2_;
    const std::vector<double>& v3r3 = this->v3rho3_;
    const std::vector<double>& v3r2s = this->v3rho2sigma_;
    const std::vector<double>& v3rs2 = this->v3rhosigma2_;
    const std::vector<double>& v3s3 = this->v3sigma3_;

    if (nspin == 1)
    {
        // Single component everywhere; all the weight sums collapse to 1.
        //
        // ... and then the whole set is scaled by 4 to reach the *singlet* normalization, the same
        // one `nspin=2` produces below and the only one the rest of the gradient code knows about.
        // libxc's unpolarized derivatives are taken w.r.t. the TOTAL density, so for a closed shell
        // (rho_u = rho_d = rho/2) an n-th derivative is 2^(n-1) smaller than the singlet spin
        // combination: d^2E/drho^2 = (f_uu+f_ud)/2 -- which is why `SpinType::S1` carries a 2 --
        // and d^3E/drho^3 = (g_uuu + 3g_uud)/4, while the nspin=2 branch below builds
        // a_s2 = g_uuu + 2g_uud + g_udd = g_uuu + 3g_uud. Hence 4 here, 2 there.
        //
        // Measured on H2/SZ/LDA: without it the GXC DMTRANS force was 0.08549 eV/Ang against the
        // nspin=2 singlet's 0.17097 -- a factor 2, being 1/4 from this and 2 from the `dm_gs`
        // channel convention (`gs_dm_channel_factor` in `lr_force.cpp`).
        constexpr double to_singlet = 4.;
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096)
#endif
        for (int i = 0;i < nrxx;++i)
        {
            dst.a_s2[i] = to_singlet * v3r3[i];
            if (!is_gga) { continue; }
            dst.a_st[i] = to_singlet * v3r2s[i] * 4.;
            dst.a_t2[i] = to_singlet * v3rs2[i] * 4.;
            dst.a_q[i] = to_singlet * v2rs[i] * 2.;
            dst.c_s[i] = to_singlet * v2rs[i] * 4.;
            dst.c_t[i] = to_singlet * v2s2[i] * 8.;
            dst.e_s2[i] = to_singlet * v3r2s[i] * 2.;
            dst.e_st[i] = to_singlet * v3rs2[i] * 8.;
            dst.e_t2[i] = to_singlet * v3s3[i] * 8.;
            dst.e_q[i] = to_singlet * v2s2[i] * 4.;
        }
        return;
    }

    // nspin=2, close shell. Every sum below runs over ORDERED spin indices, while libxc only
    // stores the unique combinations -- that is where the multiplicities come from.
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096)
#endif
    for (int i = 0;i < nrxx;++i)
    {
        const int o4 = i * 4, o6 = i * 6, o9 = i * 9, o10 = i * 10, o12 = i * 12;

        // $a_{s^2}=\sum_{\sigma\sigma'}\eta_\sigma\eta_{\sigma'}g^{\rho_u\rho_\sigma\rho_{\sigma'}}$
        // v3rho3 = (uuu, uud, udd, ddd); with the first index pinned to u the component index is
        // just the number of d's among (sigma, sigma').
        dst.a_s2[i] = eta[0] * eta[0] * v3r3[o4]
            + 2. * eta[0] * eta[1] * v3r3[o4 + 1]
            + eta[1] * eta[1] * v3r3[o4 + 2];
        if (!is_gga) { continue; }

        // ---- the local part $A$ ----
        // $a_{st}=4\sum_\sigma\eta_\sigma\sum_{\alpha\beta}\mu_{\alpha\beta}g^{\rho_u\rho_\sigma\sigma_{\alpha\beta}}$
        // v3rho2sigma = [rho-pair uu,ud,dd] x [sigma uu,ud,dd]; first rho index pinned to u means
        // rho-pair block = (sigma==d).
        double a_st = 0.;
        for (int sg = 0;sg < 2;++sg) {
            for (int b = 0;b < 3;++b) { a_st += eta[sg] * mu[b] * v3r2s[o9 + 3 * sg + b]; } }
        dst.a_st[i] = a_st * 4.;

        // $a_{t^2}=4\sum_{\alpha\beta,\gamma\delta}\mu\mu\,g^{\rho_u\sigma_{\alpha\beta}\sigma_{\gamma\delta}}$
        double a_t2 = 0.;
        for (int b = 0;b < 3;++b) {
            for (int c = 0;c < 3;++c) { a_t2 += mu[b] * mu[c] * v3rs2[o12 + p2[b][c]]; } }
        dst.a_t2[i] = a_t2 * 4.;

        // $a_q=2F$, $F=\sum_{\alpha\beta}\nu_{\alpha\beta}f^{\rho_u\sigma_{\alpha\beta}}$
        double F = 0.;
        for (int b = 0;b < 3;++b) { F += nu[b] * v2rs[o6 + b]; }
        dst.a_q[i] = F * 2.;

        // ---- the divergence part $\boldsymbol{E}$ ----
        // $P=\sum_{\alpha\beta}\theta_{\alpha\beta}\sum_{\sigma\sigma'}\eta\eta\,g^{\rho_\sigma\rho_{\sigma'}\sigma_{\alpha\beta}}$
        // rho-pair block index = number of d's in (sigma, sigma'), with multiplicity 2 for ud.
        double P = 0., Q = 0., R = 0., S = 0., T = 0., St = 0.;
        for (int a = 0;a < 3;++a)
        {
            const double w = th[a], wt = tht[a];
            if (w != 0.)
            {
                P += w * (eta[0] * eta[0] * v3r2s[o9 + 0 + a]
                    + 2. * eta[0] * eta[1] * v3r2s[o9 + 3 + a]
                    + eta[1] * eta[1] * v3r2s[o9 + 6 + a]);
                for (int sg = 0;sg < 2;++sg) {
                    for (int c = 0;c < 3;++c) { Q += w * eta[sg] * mu[c] * v3rs2[o12 + 6 * sg + p2[a][c]]; } }
                for (int c = 0;c < 3;++c) {
                    for (int d = 0;d < 3;++d) { R += w * mu[c] * mu[d] * v3s3[o10 + p3[a][c][d]]; } }
                for (int c = 0;c < 3;++c) { S += w * nu[c] * v2s2[o6 + p2[a][c]]; }
            }
            if (wt != 0.)
            {
                for (int sg = 0;sg < 2;++sg) { T += wt * eta[sg] * v2rs[o6 + 3 * sg + a]; }
                for (int c = 0;c < 3;++c) { St += wt * mu[c] * v2s2[o6 + p2[a][c]]; }
            }
        }
        dst.c_s[i] = T * 2.;
        dst.c_t[i] = St * 4.;
        dst.e_s2[i] = P;
        dst.e_st[i] = Q * 4.;
        dst.e_t2[i] = R * 4.;
        dst.e_q[i] = S * 2.;
    }
}

#endif
