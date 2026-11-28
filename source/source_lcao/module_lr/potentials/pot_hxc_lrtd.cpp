#include "pot_hxc_lrtd.h"
#include "source_base/timer.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include <set>
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include "source_lcao/module_lr/utils/lr_util.h"
#define FXC_PARA_TYPE const double* const rho, ModuleBase::matrix& v_eff, const std::vector<int>& ispin_op = { 0,0 }
namespace LR
{
    using Vec3 = ModuleBase::Vector3<double>;
    PotHxcLR::PotHxcLR(std::shared_ptr<const KernelXC> kernel,
        const ModulePW::PW_Basis& rho_basis, const UnitCell& ucell, const int nrxx,
        const int nspin, const SpinType& st)
        : PotLRBase(rho_basis, nspin, nrxx, ucell.tpiba),
          xc_kernel_components_(kernel), spin_type_(st)
    {
        if (!kernel || kernel->nspin() != nspin)
        {
            throw std::invalid_argument("PotHxcLR requires a non-null kernel with matching nspin");
        }
        this->xc_type_ = static_cast<XCType>(kernel->functional_type());
        const bool local_xc = LR_Util::has_local_xc(kernel->name());
        if (local_xc)
        {
            this->set_integral_func(this->spin_type_, this->xc_type_);
        }
        const bool gga = (this->xc_type_ == XCType::GGA || this->xc_type_ == XCType::HYB_GGA);
        this->build_spin_combos(gga);
    }

    PotHxcLR::Scratch& PotHxcLR::scratch()
    {
        static Scratch s;   // see the comment on `Scratch` in the header
        return s;
    }

    void PotHxcLR::Scratch::alloc(const int nrxx, const int npw, const int nmaxgr, const bool gga)
    {
        const std::size_t real_size = nrxx;
        const std::size_t reciprocal_size = npw;
        const std::size_t fft_size = nmaxgr;
        // resize() on an already-large vector is a no-op, so only the first call allocates.
        if (this->vxc.size() < real_size)
        {
            this->vxc.resize(real_size);
        }
        if (this->rhog.size() < reciprocal_size)
        {
            this->rhog.resize(reciprocal_size);
        }
        if (this->vg.size() < fft_size)
        {
            this->vg.resize(fft_size);
        }
        if (gga)
        {
            if (this->drho.size() < real_size)
            {
                this->drho.resize(real_size);
            }
            if (this->gdot.size() < real_size)
            {
                this->gdot.resize(real_size);
            }
        }
    }

    void PotHxcLR::build_spin_combos(const bool gga)
    {
        if (this->nspin != 2)
        {
            return;
        } // nspin=1 kernels have a single component already
        double sign = 0.;
        switch (this->spin_type_)
        {
        case SpinType::S2_singlet:
        case SpinType::S2_gs:
            sign = 1.;
            break;
        case SpinType::S2_triplet:
            sign = -1.;
            break;
        default: return;   // S2_updown indexes single components, nothing to pre-contract
        }
        const LR::KernelXC& fxc = *this->xc_kernel_components_;
        if (this->v2rho2_comb_.empty() && !fxc.v2rho2.empty())
        {
            this->v2rho2_comb_.resize(nrxx);
            const double* const f = fxc.v2rho2.data();
            double* const c = this->v2rho2_comb_.data();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (int ir = 0; ir < nrxx; ++ir)
            {
                c[ir] = f[3 * ir] + sign * f[3 * ir + 1];
            }
        }
        if (gga && this->vsigma_comb_.empty() && !fxc.vsigma.empty())
        {
            this->vsigma_comb_.resize(nrxx);
            const double* const f = fxc.vsigma.data();
            double* const c = this->vsigma_comb_.data();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (int ir = 0; ir < nrxx; ++ir)
            {
                c[ir] = 2. * f[3 * ir] + sign * f[3 * ir + 1];
            }
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
        for (int ir = 0; ir < nrxx; ++ir)
        {
            v[ir] += factor * vg[ir].real();
        }
    }

    void PotHxcLR::cal_v_eff(double** rho, const UnitCell& ucell, ModuleBase::matrix& v_eff, const std::vector<int>& ispin_op) const
    {
        ModuleBase::TITLE("PotHxcLR", "cal_v_eff");
        ModuleBase::timer::start("PotHxcLR", "cal_v_eff");

        // Hartree prefactor (0 = this spin combination has no Hartree term at all)
        double hartree_factor = 0.;
        switch (this->spin_type_)
        {
        case SpinType::S1_gs:
        case SpinType::S2_updown:
        case SpinType::S2_gs:
            hartree_factor = 1.;
            break;
        case SpinType::S1:
        case SpinType::S2_singlet:
            hartree_factor = 2.;
            break;
        default: break;
        }
        const bool local_xc = LR_Util::has_local_xc(this->xc_kernel_components_->name());
        const bool gga = local_xc && (this->xc_type_ == XCType::GGA || this->xc_type_ == XCType::HYB_GGA);

        // $\rho^X(G)$ is needed by the Hartree term and by the GGA gradient, and used to be
        // transformed twice (once inside `H_Hartree_pw::v_hartree`, once inside `LR_Util::grad`).
        if (hartree_factor != 0. || gga)
        {
            scratch().alloc(nrxx, this->rho_basis_.npw, this->rho_basis_.nmaxgr, gga);
            this->rho_basis_.real2recip(rho[0], scratch().rhog.data());
        }
        if (hartree_factor != 0.)
        {
            this->add_v_hartree(ucell, v_eff, hartree_factor);
        }

        // XC
        if (!local_xc)
        {
            ModuleBase::timer::end("PotHxcLR", "cal_v_eff");
            return;
        }
#ifdef __LIBXC
        this->kernel_to_potential_.at(spin_type_)(rho[0], v_eff, ispin_op);
#else
        throw std::domain_error("PotHxcLR XC type =" + std::to_string(this->xc_type_)
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
        std::map<SpinType, Tfunc>& funcs = this->kernel_to_potential_;
        const LR::KernelXC& fxc = *this->xc_kernel_components_;
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
                        for (int ir = 0; ir < nrxx; ++ir)
                        {
                            v[ir] += fac * f[ir] * rho[ir];
                        }
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
                        for (int ir = 0; ir < nrxx; ++ir)
                        {
                            v[ir] += fac * f[ir] * rho[ir];
                        }
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
                        for (int ir = 0; ir < nrxx; ++ir)
                        {
                            v[ir] += ModuleBase::e2 * f[3 * ir + is] * rho[ir];
                        }
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
                        const auto& coefficients = fxc.closed_shell_gga(false);
                        const double* const rs2 = coefficients.rho_sigma.data();
                        const double* const ss4 = coefficients.sigma_sigma.data();
                        const Vec3* const dgs = fxc.drho_gs.at(0).data();
                        const double* const vsig = fxc.vsigma.data();
                        const double* const rr = fxc.v2rho2.data();

                        //1. $\partial E/\partial\rho = 2f^{\rho\sigma}*\nabla\rho*\rho_1+4f^{\sigma\sigma}\nabla\rho(\nabla\rho\cdot\nabla\rho_1)+2v^\sigma\nabla\rho_1$
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
                        for (int ir = 0;ir < nrxx;++ir)
                        {
                            e_drho[ir] = -(dgs[ir] * (rs2[ir] * rho[ir]
                                + ss4[ir] * (dgs[ir] * drho[ir]))
                                + drho[ir] * (vsig[ir] * 2.));
                        }
                        XC_Functional::grad_dot(e_drho, vxc, &this->rho_basis_, this->tpiba_);

                        // 2. $f^{\rho\rho}\rho_1+2f^{\rho\sigma}\nabla\rho\cdot\nabla\rho_1$
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
                        for (int ir = 0;ir < nrxx;++ir)
                        {
                            vxc[ir] += rr[ir] * rho[ir] + rs2[ir] * (dgs[ir] * drho[ir]);
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
                        const auto& coefficients = fxc.closed_shell_gga(triplet);
                        const double* const rs = coefficients.rho_sigma.data();
                        const double* const ss = coefficients.sigma_sigma.data();
                        const Vec3* const dgs = fxc.drho_gs.at(0).data();
                        const double* const vsig = this->vsigma_comb_.data();   // 2*vsigma_uu -+ vsigma_ud
                        const double* const rr = this->v2rho2_comb_.data();     // v2rho2_uu -+ v2rho2_ud

                        // 1. the terms in grad_dot int f_uu
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
                        for (int ir = 0;ir < nrxx;++ir)
                        {
                            gdot[ir] = -(dgs[ir] * (rs[ir] * rho[ir]
                                + ss[ir] * (dgs[ir] * drho[ir]))
                                + drho[ir] * vsig[ir]);
                        }
                        XC_Functional::grad_dot(gdot, vxc, &this->rho_basis_, this->tpiba_);

                        // 2. terms not in grad_dot
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
                        for (int ir = 0;ir < nrxx;++ir)
                        {
                            vxc[ir] += rho[ir] * rr[ir] + (drho[ir] * dgs[ir]) * rs[ir];
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

                        const int output_spin = ispin_op[0];
                        const int input_spin = ispin_op[1];
                        if (output_spin < 0 || output_spin > 1 || input_spin < 0 || input_spin > 1)
                        {
                            throw std::runtime_error("Invalid ispin_op");
                        }
                        const double* const rho_sigma = fxc.v2rhosigma.data();
                        const double* const sigma_sigma = fxc.v2sigma2.data();
                        // Each input spin perturbs two adjacent sigma invariants. Select all
                        // spin-dependent gradients and coefficient indices before the grid loops.
                        const int in0 = input_spin;
                        const int in1 = input_spin + 1;
                        const int out0 = output_spin;
                        const int out1 = output_spin + 1;
                        const Vec3* const gs_gradient[2] = {dgs_u, dgs_d};
                        const Vec3* const in_gradient0 = gs_gradient[libxc_idx::chan[input_spin][in0]];
                        const Vec3* const in_gradient1 = gs_gradient[libxc_idx::chan[input_spin][in1]];
                        const Vec3* const out_gradient0 = gs_gradient[libxc_idx::chan[output_spin][out0]];
                        const Vec3* const out_gradient1 = gs_gradient[libxc_idx::chan[output_spin][out1]];
                        const double in_weight0 = libxc_idx::theta[input_spin][in0];
                        const double in_weight1 = libxc_idx::theta[input_spin][in1];
                        const double out_weight0 = libxc_idx::theta[output_spin][out0];
                        const double out_weight1 = libxc_idx::theta[output_spin][out1];
                        const int rs0 = 3 * input_spin + out0;
                        const int rs1 = 3 * input_spin + out1;
                        const int ss00 = libxc_idx::p2[out0][in0];
                        const int ss01 = libxc_idx::p2[out0][in1];
                        const int ss10 = libxc_idx::p2[out1][in0];
                        const int ss11 = libxc_idx::p2[out1][in1];
                        const int local_rs0 = 3 * output_spin + in0;
                        const int local_rs1 = 3 * output_spin + in1;
                        const int spin_pair = output_spin + input_spin;
                        const double gradient_weight = (output_spin == input_spin) ? 2. : 1.;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
                        for (int ir = 0; ir < nrxx; ++ir)
                        {
                            const double sigma0 = in_weight0 * (in_gradient0[ir] * drho[ir]);
                            const double sigma1 = in_weight1 * (in_gradient1[ir] * drho[ir]);
                            const double* const rs = rho_sigma + 6 * ir;
                            const double* const ss = sigma_sigma + 6 * ir;
                            const double derivative0 = rs[rs0] * rho[ir] + ss[ss00] * sigma0 + ss[ss01] * sigma1;
                            const double derivative1 = rs[rs1] * rho[ir] + ss[ss10] * sigma0 + ss[ss11] * sigma1;
                            gdot[ir] = -(out_gradient0[ir] * (out_weight0 * derivative0)
                                + out_gradient1[ir] * (out_weight1 * derivative1)
                                + drho[ir] * (gradient_weight * vsig[3 * ir + spin_pair]));
                        }
                        XC_Functional::grad_dot(gdot, vxc, &this->rho_basis_, this->tpiba_);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
                        for (int ir = 0; ir < nrxx; ++ir)
                        {
                            const double sigma0 = in_weight0 * (in_gradient0[ir] * drho[ir]);
                            const double sigma1 = in_weight1 * (in_gradient1[ir] * drho[ir]);
                            const double* const rs = rho_sigma + 6 * ir;
                            vxc[ir] += rho[ir] * rr[3 * ir + spin_pair]
                                + rs[local_rs0] * sigma0 + rs[local_rs1] * sigma1;
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
            throw std::domain_error("PotHxcLR XC type =" + std::to_string(this->xc_type_)
                + " unfinished in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
        }
    }
} // namespace LR
