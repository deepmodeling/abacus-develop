#include "pot_hxc_lrtd.h"
#include "source_io/module_parameter/parameter.h"
#include "source_estate/module_pot/h_hartree_pw.h"
#include "source_base/timer.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include <set>
#include <algorithm>
#include "source_lcao/module_lr/utils/lr_util.h"
#include "source_lcao/module_lr/utils/lr_util_xc.hpp"
#define FXC_PARA_TYPE const double* const rho, ModuleBase::matrix& v_eff, const std::vector<int>& ispin_op = { 0,0 }
namespace LR
{
    using Vec3 = ModuleBase::Vector3<double>;
    /// the `nspin` `KernelXC` is built with: 1 for a non-magnetic calculation, 2 otherwise.
    /// Kept next to `PotLRBase`'s own expression so `make_kernel` cannot drift away from it.
    static int kernel_nspin()
    {
        return (PARAM.inp.nspin == 1 || (PARAM.inp.nspin == 4 && !PARAM.globalv.domag && !PARAM.globalv.domag_z)) ? 1 : 2;
    }

    std::shared_ptr<const KernelXC> PotHxcLR::make_kernel(const std::string& xc_kernel,
        const ModulePW::PW_Basis& rho_basis, const UnitCell& ucell, const Charge& chg_gs,
        const Parallel_Grid& pgrid, const bool openshell, const int gxc_spin,
        const std::vector<std::string>& lr_init_xc_kernel)
    {   //calls XC_Functional::set_func_type and libxc
        return std::make_shared<const KernelXC>(rho_basis, ucell, chg_gs, pgrid, kernel_nspin(),
            xc_kernel, lr_init_xc_kernel, openshell, gxc_spin);
    }

    // constructor for exchange-correlation kernel
    PotHxcLR::PotHxcLR(const std::string& xc_kernel, const ModulePW::PW_Basis& rho_basis, const UnitCell& ucell,
        const Charge& chg_gs/*ground state*/, const Parallel_Grid& pgrid,
        const SpinType& st, const std::vector<std::string>& lr_init_xc_kernel, const int gxc_spin)
        :PotLRBase(rho_basis, kernel_nspin(), chg_gs.nrxx, ucell.tpiba),
        xc_kernel_(xc_kernel), spin_type_(st),
        pot_hartree_(LR_Util::make_unique<elecstate::PotHartree>(&rho_basis)),
        xc_kernel_components_(make_kernel(xc_kernel, rho_basis, ucell, chg_gs, pgrid, (st == SpinType::S2_updown), gxc_spin, lr_init_xc_kernel)),
        xc_type_(XCType(XC_Functional::get_func_type()))
    {
        if (LR_Util::has_local_xc(xc_kernel)) { this->set_integral_func(this->spin_type_, this->xc_type_); }
    }

    PotHxcLR::PotHxcLR(std::shared_ptr<const KernelXC> kernel, const std::string& xc_kernel,
        const ModulePW::PW_Basis& rho_basis, const UnitCell& ucell, const int nrxx, const SpinType& st)
        :PotLRBase(rho_basis, kernel_nspin(), nrxx, ucell.tpiba),
        xc_kernel_(xc_kernel), spin_type_(st),
        pot_hartree_(LR_Util::make_unique<elecstate::PotHartree>(&rho_basis)),
        xc_kernel_components_(std::move(kernel)),
        xc_type_(XCType(XC_Functional::get_func_type()))
    {
        assert(this->xc_kernel_components_ != nullptr);
        if (LR_Util::has_local_xc(xc_kernel)) { this->set_integral_func(this->spin_type_, this->xc_type_); }
    }


    PotHxcLR::Scratch& PotHxcLR::scratch()
    {
        static Scratch s;   // see the comment on `Scratch` in the header
        return s;
    }

    void PotHxcLR::Scratch::alloc(const int nrxx, const int npw, const int nmaxgr, const bool gga)
    {
        // resize() on an already-large vector is a no-op, so only the first call allocates.
        if (static_cast<int>(this->vxc.size()) < nrxx) { this->vxc.resize(nrxx); }
        if (static_cast<int>(this->rhog.size()) < npw) { this->rhog.resize(npw); }
        if (static_cast<int>(this->vg.size()) < nmaxgr) { this->vg.resize(nmaxgr); }
        if (gga)
        {
            if (static_cast<int>(this->drho.size()) < nrxx) { this->drho.resize(nrxx); }
            if (static_cast<int>(this->gdot.size()) < nrxx) { this->gdot.resize(nrxx); }
        }
    }

    void PotHxcLR::build_spin_combos(const bool gga) const
    {
        if (this->nspin != 2) { return; }   // nspin=1 kernels have a single component already
        double sign = 0.;
        switch (this->spin_type_)
        {
        case SpinType::S2_singlet: case SpinType::S2_gs: sign = 1.; break;
        case SpinType::S2_triplet:                       sign = -1.; break;
        default: return;   // S2_updown indexes single components, nothing to pre-contract
        }
        auto& fxc = *this->xc_kernel_components_;
        if (this->v2rho2_comb_.empty() && !fxc.v2rho2.empty())
        {
            this->v2rho2_comb_.resize(nrxx);
            const double* const f = fxc.v2rho2.data();
            double* const c = this->v2rho2_comb_.data();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (int ir = 0; ir < nrxx; ++ir) { c[ir] = f[3 * ir] + sign * f[3 * ir + 1]; }
        }
        if (gga && this->vsigma_comb_.empty() && !fxc.vsigma.empty())
        {
            this->vsigma_comb_.resize(nrxx);
            const double* const f = fxc.vsigma.data();
            double* const c = this->vsigma_comb_.data();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (int ir = 0; ir < nrxx; ++ir) { c[ir] = 2. * f[3 * ir] + sign * f[3 * ir + 1]; }
        }
    }

    void PotHxcLR::add_v_hartree(const UnitCell& ucell, ModuleBase::matrix& v_eff, const double factor) const
    {
        const int npw = this->rho_basis_.npw;
        const std::complex<double>* const rhog = scratch().rhog.data();
        std::complex<double>* const vg = scratch().vg.data();
        const int ig0 = this->rho_basis_.ig_gge0;
        const double pre = ModuleBase::e2 * ModuleBase::FOUR_PI / ucell.tpiba2;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ig = 0; ig < npw; ++ig)
        {
            vg[ig] = (ig == ig0) ? std::complex<double>(0., 0.)   // V(G=0) = 0
                                 : (pre / this->rho_basis_.gg[ig]) * rhog[ig];
        }
        this->rho_basis_.recip2real(vg, vg);
        double* const v = v_eff.c;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ir = 0; ir < nrxx; ++ir) { v[ir] += factor * vg[ir].real(); }
    }

    void PotHxcLR::cal_v_eff(double** rho, const UnitCell& ucell, ModuleBase::matrix& v_eff, const std::vector<int>& ispin_op) const
    {
        ModuleBase::TITLE("PotHxcLR", "cal_v_eff");
        ModuleBase::timer::start("PotHxcLR", "cal_v_eff");

        // Hartree prefactor (0 = this spin combination has no Hartree term at all)
        double hartree_factor = 0.;
        switch (this->spin_type_)
        {
        case SpinType::S1_gs: case SpinType::S2_updown: case SpinType::S2_gs: hartree_factor = 1.; break;
        case SpinType::S1: case SpinType::S2_singlet:                         hartree_factor = 2.; break;
        default: break;
        }
        const bool local_xc = LR_Util::has_local_xc(this->xc_kernel_);
        const bool gga = local_xc && (this->xc_type_ == XCType::GGA || this->xc_type_ == XCType::HYB_GGA);

        // $\rho^X(G)$ is needed by the Hartree term and by the GGA gradient, and used to be
        // transformed twice (once inside `H_Hartree_pw::v_hartree`, once inside `LR_Util::grad`).
        if (hartree_factor != 0. || gga)
        {
            scratch().alloc(nrxx, this->rho_basis_.npw, this->rho_basis_.nmaxgr, gga);
            this->rho_basis_.real2recip(rho[0], scratch().rhog.data());
        }
        if (hartree_factor != 0.) { this->add_v_hartree(ucell, v_eff, hartree_factor); }

        // XC
        if (!local_xc)
        {
            ModuleBase::timer::end("PotHxcLR", "cal_v_eff");
            return;
        }
#ifdef __LIBXC
        this->build_spin_combos(gga);
        this->kernel_to_potential_.at(spin_type_)(rho[0], v_eff, ispin_op);
#else
        throw std::domain_error("GlobalV::XC_Functional::get_func_type() =" + std::to_string(XC_Functional::get_func_type())
            + " unfinished in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
#endif
        ModuleBase::timer::end("PotHxcLR", "cal_v_eff");
    }

    // In every integrand below the kernel arrays are read through raw pointers hoisted out of
    // the loop, and the loops carry an OpenMP directive. Both matter: these are pure
    // element-wise passes over nrxx that used to run single-threaded, with `.at()` (a bounds
    // check per access, which also blocks vectorization) on every kernel component.
    void PotHxcLR::set_integral_func(const SpinType& s, const XCType& xc)
    {
        auto& funcs = this->kernel_to_potential_;
        auto& fxc = *this->xc_kernel_components_;
        if (xc == XCType::LDA) {
            switch (s)
            {
            case SpinType::S1:
            case SpinType::S1_gs:
            {
                // S1_gs is exactly half of S1 (see the SpinType doc in the header).
                const double prefac = (s == SpinType::S1_gs) ? 1.0 : 2.0;
                funcs[s] = [this, &fxc, prefac](FXC_PARA_TYPE)->void
                    {
                        const double* const f = fxc.v2rho2.data();
                        double* const v = v_eff.c;
                        const double fac = ModuleBase::e2 * prefac;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
                        for (int ir = 0;ir < nrxx;++ir) { v[ir] += fac * f[ir] * rho[ir]; }
                    };
                break;
            }
            case SpinType::S2_singlet:
            case SpinType::S2_gs:
            case SpinType::S2_triplet:
            {
                // S2_gs is exactly half of S2_singlet (see the SpinType doc in the header).
                // The singlet/triplet difference is entirely in `v2rho2_comb_`'s sign.
                const double prefac = (s == SpinType::S2_gs) ? 0.5 : 1.0;
                funcs[s] = [this, prefac](FXC_PARA_TYPE)->void
                    {
                        const double* const f = this->v2rho2_comb_.data();
                        double* const v = v_eff.c;
                        const double fac = ModuleBase::e2 * prefac;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
                        for (int ir = 0;ir < nrxx;++ir) { v[ir] += fac * f[ir] * rho[ir]; }
                    };
                break;
            }
            case SpinType::S2_updown:
                funcs[s] = [this, &fxc](FXC_PARA_TYPE)->void
                    {
                        assert(ispin_op.size() >= 2);
                        const int is = ispin_op[0] + ispin_op[1];
                        const double* const f = fxc.v2rho2.data();
                        double* const v = v_eff.c;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
                        for (int ir = 0;ir < nrxx;++ir) { v[ir] += ModuleBase::e2 * f[3 * ir + is] * rho[ir]; }
                    };
                break;
            default:
                throw std::domain_error("SpinType =" + std::to_string(static_cast<int>(s))
                    + " unfinished in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
                break;
            }
        }
        else if (xc == XCType::GGA || xc == XCType::HYB_GGA) {
            switch (s)
            {
            case SpinType::S1:
            case SpinType::S1_gs:
            {
                const double prefac = (s == SpinType::S1_gs) ? 1.0 : 2.0;   // S1_gs is exactly half of S1.
                funcs[s] = [this, &fxc, prefac](FXC_PARA_TYPE)->void
                    {
                        // transition density gradient, from the $\rho^X(G)$ `cal_v_eff` already built
                        Vec3* const drho = scratch().drho.data();
                        XC_Functional::grad_rho(scratch().rhog.data(), drho, &this->rho_basis_, this->tpiba_);

                        double* const vxc = scratch().vxc.data();
                        Vec3* const e_drho = scratch().gdot.data();
                        const Vec3* const rs2 = fxc.v2rhosigma_2drho.data();
                        const Vec3* const ss4 = fxc.v2sigma2_4drho.data();
                        const Vec3* const dgs = fxc.drho_gs.at(0).data();
                        const double* const vsig = fxc.vsigma.data();
                        const double* const rr = fxc.v2rho2.data();

                        //1. $\partial E/\partial\rho = 2f^{\rho\sigma}*\nabla\rho*\rho_1+4f^{\sigma\sigma}\nabla\rho(\nabla\rho\cdot\nabla\rho_1)+2v^\sigma\nabla\rho_1$
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
                        for (int ir = 0;ir < nrxx;++ir)
                        {
                            e_drho[ir] = -(rs2[ir] * rho[ir]
                                + ss4[ir] * (dgs[ir] * drho[ir])
                                + drho[ir] * (vsig[ir] * 2.));
                        }
                        XC_Functional::grad_dot(e_drho, vxc, &this->rho_basis_, this->tpiba_);

                        // 2. $f^{\rho\rho}\rho_1+2f^{\rho\sigma}\nabla\rho\cdot\nabla\rho_1$
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
                        for (int ir = 0;ir < nrxx;++ir)
                        {
                            vxc[ir] += rr[ir] * rho[ir] + rs2[ir] * drho[ir];
                        }
                        BlasConnector::axpy(nrxx, ModuleBase::e2 * prefac, vxc, 1, v_eff.c, 1);
                    };
                break;
            }
            case SpinType::S2_singlet:
            case SpinType::S2_gs:
            case SpinType::S2_triplet:
            {
                // S2_gs is exactly half of S2_singlet; the whole expression is linear in the
                // kernel, so scaling the final axpy is enough. Singlet and triplet differ only
                // in which pre-contracted kernel arrays are read, so they share one lambda.
                const double prefac = (s == SpinType::S2_gs) ? 0.5 : 1.0;
                const bool triplet = (s == SpinType::S2_triplet);
                funcs[s] = [this, &fxc, prefac, triplet](FXC_PARA_TYPE)->void
                    {
                        Vec3* const drho = scratch().drho.data();
                        XC_Functional::grad_rho(scratch().rhog.data(), drho, &this->rho_basis_, this->tpiba_);

                        double* const vxc = scratch().vxc.data();
                        Vec3* const gdot = scratch().gdot.data();
                        const Vec3* const rs = triplet ? fxc.v2rhosigma_drho_triplet.data()
                                                       : fxc.v2rhosigma_drho_singlet.data();
                        const Vec3* const ss = triplet ? fxc.v2sigma2_drho_triplet.data()
                                                       : fxc.v2sigma2_drho_singlet.data();
                        const Vec3* const dgs = fxc.drho_gs.at(0).data();
                        const double* const vsig = this->vsigma_comb_.data();   // 2*vsigma_uu -+ vsigma_ud
                        const double* const rr = this->v2rho2_comb_.data();     // v2rho2_uu -+ v2rho2_ud

                        // 1. the terms in grad_dot int f_uu
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
                        for (int ir = 0;ir < nrxx;++ir)
                        {
                            gdot[ir] = -(rs[ir] * rho[ir]
                                + ss[ir] * (dgs[ir] * drho[ir])
                                + drho[ir] * vsig[ir]);
                        }
                        XC_Functional::grad_dot(gdot, vxc, &this->rho_basis_, this->tpiba_);

                        // 2. terms not in grad_dot
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
                        for (int ir = 0;ir < nrxx;++ir)
                        {
                            vxc[ir] += rho[ir] * rr[ir] + drho[ir] * rs[ir];
                        }
                        BlasConnector::axpy(nrxx, ModuleBase::e2 * prefac, vxc, 1, v_eff.c, 1);
                    };
                break;
            }
            case SpinType::S2_updown:
                funcs[s] = [this, &fxc](FXC_PARA_TYPE)->void
                    {
                        assert(ispin_op.size() >= 2);
                        Vec3* const drho = scratch().drho.data();
                        XC_Functional::grad_rho(scratch().rhog.data(), drho, &this->rho_basis_, this->tpiba_);

                        double* const vxc = scratch().vxc.data();
                        Vec3* const gdot = scratch().gdot.data();
                        const Vec3* const dgs_u = fxc.drho_gs.at(0).data();
                        const Vec3* const dgs_d = fxc.drho_gs.at(1).data();
                        const double* const vsig = fxc.vsigma.data();
                        const double* const rr = fxc.v2rho2.data();

                        // the four (sigma, sigma') blocks: which pre-contracted arrays to read,
                        // which vsigma component multiplies drho^X, and which v2rho2 component.
                        const int blk = ispin_op[0] << 1 | ispin_op[1];
                        const Vec3* rs = nullptr; const Vec3* su = nullptr; const Vec3* sd = nullptr;
                        int i_vsig = 0, i_rr = 0; double vsig_fac = 1.;
                        // `rs2` is the v2rhosigma array used again in the non-divergence term;
                        // for the off-diagonal blocks it is the transposed one (ud vs du).
                        const Vec3* rs2 = nullptr;
                        switch (blk)
                        {
                        case 0: // (0,0)
                            rs = fxc.v2rhosigma_drho_uu.data(); rs2 = rs;
                            su = fxc.v2sigma2_drho_uu_u.data(); sd = fxc.v2sigma2_drho_uu_d.data();
                            i_vsig = 0; vsig_fac = 2.; i_rr = 0;
                            break;
                        case 1: // (0,1): rho_d, drho_u
                            rs = fxc.v2rhosigma_drho_du.data(); rs2 = fxc.v2rhosigma_drho_ud.data();
                            su = fxc.v2sigma2_drho_ud_u.data(); sd = fxc.v2sigma2_drho_ud_d.data();
                            i_vsig = 1; vsig_fac = 1.; i_rr = 1;
                            break;
                        case 2: // (1,0): rho_u, drho_d
                            rs = fxc.v2rhosigma_drho_ud.data(); rs2 = fxc.v2rhosigma_drho_du.data();
                            su = fxc.v2sigma2_drho_du_u.data(); sd = fxc.v2sigma2_drho_du_d.data();
                            i_vsig = 1; vsig_fac = 1.; i_rr = 1;
                            break;
                        case 3: // (1,1)
                            rs = fxc.v2rhosigma_drho_dd.data(); rs2 = rs;
                            su = fxc.v2sigma2_drho_dd_u.data(); sd = fxc.v2sigma2_drho_dd_d.data();
                            i_vsig = 2; vsig_fac = 2.; i_rr = 2;
                            break;
                        default:
                            throw std::runtime_error("Invalid ispin_op");
                        }
                        // 1. the terms in grad_dot int f_uu
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
                        for (int ir = 0;ir < nrxx;++ir)
                        {
                            gdot[ir] = -(rs[ir] * rho[ir]
                                + su[ir] * (dgs_u[ir] * drho[ir])
                                + sd[ir] * (dgs_d[ir] * drho[ir])
                                + drho[ir] * (vsig[ir * 3 + i_vsig] * vsig_fac));
                        }
                        XC_Functional::grad_dot(gdot, vxc, &this->rho_basis_, this->tpiba_);

                        // 2. terms not in grad_dot
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
                        for (int ir = 0;ir < nrxx;++ir)
                        {
                            vxc[ir] += rho[ir] * rr[3 * ir + i_rr] + drho[ir] * rs2[ir];
                        }
                        BlasConnector::axpy(nrxx, ModuleBase::e2, vxc, 1, v_eff.c, 1);
                    };
                break;
            default:
                throw std::domain_error("SpinType =" + std::to_string(static_cast<int>(s)) + "for GGA or HYB_GGA is unfinished in "
                    + std::string(__FILE__) + " line " + std::to_string(__LINE__));
                break;
            }
        }
        else
        {
            throw std::domain_error("GlobalV::XC_Functional::get_func_type() =" + std::to_string(XC_Functional::get_func_type())
                + " unfinished in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
        }
    }
} // namespace LR
