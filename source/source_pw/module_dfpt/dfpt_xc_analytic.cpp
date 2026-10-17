#ifdef __LIBXC

#include "dfpt_xc_analytic.h"

#include "source_base/constants.h"
#include "source_base/tool_quit.h"

#include <xc.h>

#include <cmath>
#include <complex>
#include <vector>

namespace ModuleDFPT
{

namespace
{
// QE dgradcorr thresholds: no gradient-kernel response where the
// ground-state density or its gradient is negligible
const double epsr = 1.0e-6;
const double epsg = 1.0e-10;

// libxc works in Hartree; the DFPT potentials are in Ry (e2 = 2)
const double ha_to_ry = 2.0;
} // namespace

XC_First_Order_Analytic::XC_First_Order_Analytic(const std::vector<int>& func_id,
                                                 const int func_type,
                                                 const std::vector<double>& rho0_r,
                                                 ModulePW::PW_Basis* pw_rho)
    : func_type_(func_type), func_id_(func_id), pw_rho_(pw_rho)
{
    ModuleBase::TITLE("XC_First_Order_Analytic", "XC_First_Order_Analytic");
    if (func_type != 1 && func_type != 2)
    {
        ModuleBase::WARNING_QUIT("XC_First_Order_Analytic", "only LDA (1) and GGA (2) kernels are supported");
    }
    if (func_id_.empty() || pw_rho_ == nullptr || static_cast<int>(rho0_r.size()) != pw_rho_->nrxx)
    {
        ModuleBase::WARNING_QUIT("XC_First_Order_Analytic", "invalid func_id / basis / rho0 size");
    }
    this->init_gradient_(rho0_r);
    this->init_kernel_fields_(rho0_r);
}

void XC_First_Order_Analytic::set_context(const ModuleBase::Vector3<double>& q_cart) const
{
    q_cart_ = q_cart;
    q_set_ = true;
}

void XC_First_Order_Analytic::init_gradient_(const std::vector<double>& rho0_r)
{
    const int nrxx = pw_rho_->nrxx;
    const int npw = pw_rho_->npw;
    std::vector<std::complex<double>> rho0_g(npw);
    pw_rho_->real2recip(rho0_r.data(), rho0_g.data());
    grho_.assign(3, std::vector<double>(nrxx, 0.0));
    std::vector<std::complex<double>> g_g(npw);
    for (int ipol = 0; ipol < 3; ++ipol)
    {
        for (int ig = 0; ig < npw; ++ig)
        {
            // pw gcar is dimensionless (2 pi not folded in): the physical
            // wavevector is tpiba * gcar, in 1/lat0
            const double g = pw_rho_->tpiba * pw_rho_->gcar[ig][ipol];
            g_g[ig] = std::complex<double>(0.0, g) * rho0_g[ig];
        }
        std::vector<std::complex<double>> g_r(nrxx);
        pw_rho_->recip2real(g_g.data(), g_r.data());
        for (int ir = 0; ir < nrxx; ++ir)
        {
            grho_[ipol][ir] = g_r[ir].real();
        }
    }
}

void XC_First_Order_Analytic::init_kernel_fields_(const std::vector<double>& rho0_r)
{
    const int nrxx = pw_rho_->nrxx;
    k_rr_.assign(nrxx, 0.0);
    k_sr_.assign(nrxx, 0.0);
    k_ss_.assign(nrxx, 0.0);
    k_s_.assign(nrxx, 0.0);

    // sigma = |grad rho0|^2 for the GGA branch
    std::vector<double> sigma(nrxx, 0.0);
    if (func_type_ == 2)
    {
        for (int ir = 0; ir < nrxx; ++ir)
        {
            for (int ipol = 0; ipol < 3; ++ipol)
            {
                sigma[ir] += grho_[ipol][ir] * grho_[ipol][ir];
            }
        }
    }

    // open every functional component once, accumulate the fields, close
    for (size_t ifunc = 0; ifunc < func_id_.size(); ++ifunc)
    {
        xc_func_type func;
        if (xc_func_init(&func, func_id_[ifunc], XC_UNPOLARIZED) != 0)
        {
            ModuleBase::WARNING_QUIT("XC_First_Order_Analytic", "xc_func_init failed");
        }
        if (func_type_ == 1)
        {
            for (int ir = 0; ir < nrxx; ++ir)
            {
                double v2rho2 = 0.0;
                const double rho = rho0_r[ir];
                xc_lda_fxc(&func, 1, &rho, &v2rho2);
                // libxc vrho is already the potential d(rho e)/d rho, so
                // the kernel is its plain derivative (QE dmuxc convention)
                k_rr_[ir] += ha_to_ry * v2rho2;
            }
        }
        else
        {
            for (int ir = 0; ir < nrxx; ++ir)
            {
                const double rho = rho0_r[ir];
                const double sig = sigma[ir];
                if (std::abs(rho) <= epsr || sig <= epsg)
                {
                    continue; // QE dgradcorr guard: no kernel response here
                }
                double zk = 0.0;
                double vrho = 0.0;
                double vsigma = 0.0;
                xc_gga_exc_vxc(&func, 1, &rho, &sig, &zk, &vrho, &vsigma);
                double v2rho2 = 0.0;
                double v2rhosigma = 0.0;
                double v2sigma2 = 0.0;
                xc_gga_fxc(&func, 1, &rho, &sig, &v2rho2, &v2rhosigma, &v2sigma2);
                // sigma = |grad rho0|^2 chain factors of the dgradcorr
                // algebra: s1 = grho.gdrho = d sigma / 2, so the response
                // e_rho_sigma d(sigma) carries the factor 2 and the
                // e_sigma_sigma term the factor 4 (QE xc_wrapper_d_gga
                // convention; no rho renormalization, the potential-form
                // kernel cancels it)
                k_rr_[ir] += ha_to_ry * v2rho2;
                k_sr_[ir] += ha_to_ry * (2.0 * v2rhosigma);
                k_ss_[ir] += ha_to_ry * (4.0 * v2sigma2);
                k_s_[ir] += ha_to_ry * (2.0 * vsigma);
            }
        }
        xc_func_end(&func);
    }
}

void XC_First_Order_Analytic::apply(const std::vector<std::complex<double>>& drho_r,
                                    std::vector<std::complex<double>>& dvxc_r) const
{
    const int nrxx = pw_rho_->nrxx;
    const int npw = pw_rho_->npw;
    if (static_cast<int>(drho_r.size()) != nrxx)
    {
        ModuleBase::WARNING_QUIT("XC_First_Order_Analytic", "drho_r is not on the rho grid");
    }
    dvxc_r.assign(nrxx, std::complex<double>(0.0, 0.0));
    if (func_type_ == 1)
    {
        // LDA: dvxc = f_rho A, no q dependence
        for (int ir = 0; ir < nrxx; ++ir)
        {
            dvxc_r[ir] = k_rr_[ir] * drho_r[ir];
        }
        return;
    }

    // GGA: dgradcorr with the total-field gradient i (G + q)
    const double qx = q_set_ ? q_cart_.x : 0.0;
    const double qy = q_set_ ? q_cart_.y : 0.0;
    const double qz = q_set_ ? q_cart_.z : 0.0;
    std::vector<std::complex<double>> a_g(npw);
    pw_rho_->real2recip(drho_r.data(), a_g.data());

    // gdrho(ipol, r) = F^-1[ i tpiba (G + q)_ipol A(G) ]: physical
    // gradient of the carrier-carrying field, 1/lat0 units
    std::vector<std::vector<std::complex<double>>> gdrho(3, std::vector<std::complex<double>>(nrxx, std::complex<double>(0.0, 0.0)));
    std::vector<std::complex<double>> g_g(npw);
    for (int ipol = 0; ipol < 3; ++ipol)
    {
        const double q_ipol = ipol == 0 ? qx : (ipol == 1 ? qy : qz);
        for (int ig = 0; ig < npw; ++ig)
        {
            const double w = pw_rho_->tpiba * (pw_rho_->gcar[ig][ipol] + q_ipol);
            g_g[ig] = std::complex<double>(0.0, w) * a_g[ig];
        }
        pw_rho_->recip2real(g_g.data(), gdrho[ipol].data());
    }

    // s1 = grad rho0 . gdrho
    std::vector<std::complex<double>> s1(nrxx, std::complex<double>(0.0, 0.0));
    for (int ir = 0; ir < nrxx; ++ir)
    {
        for (int ipol = 0; ipol < 3; ++ipol)
        {
            s1[ir] += grho_[ipol][ir] * gdrho[ipol][ir];
        }
    }

    // local part and the h field whose q-gradient completes the response
    std::vector<std::vector<std::complex<double>>> h(3, std::vector<std::complex<double>>(nrxx, std::complex<double>(0.0, 0.0)));
    for (int ir = 0; ir < nrxx; ++ir)
    {
        dvxc_r[ir] = k_rr_[ir] * drho_r[ir] + k_sr_[ir] * s1[ir];
        const std::complex<double> fac = k_sr_[ir] * drho_r[ir] + k_ss_[ir] * s1[ir];
        for (int ipol = 0; ipol < 3; ++ipol)
        {
            h[ipol][ir] = fac * grho_[ipol][ir] + k_s_[ir] * gdrho[ipol][ir];
        }
    }

    // dvxc -= F^-1[ i tpiba (G + q) . h(G) ]
    std::vector<std::complex<double>> dh_g(npw, std::complex<double>(0.0, 0.0));
    std::vector<std::complex<double>> h_g(npw);
    for (int ipol = 0; ipol < 3; ++ipol)
    {
        const double q_ipol = ipol == 0 ? qx : (ipol == 1 ? qy : qz);
        pw_rho_->real2recip(h[ipol].data(), h_g.data());
        for (int ig = 0; ig < npw; ++ig)
        {
            const double w = pw_rho_->tpiba * (pw_rho_->gcar[ig][ipol] + q_ipol);
            dh_g[ig] += std::complex<double>(0.0, w) * h_g[ig];
        }
    }
    std::vector<std::complex<double>> dh_r(nrxx);
    pw_rho_->recip2real(dh_g.data(), dh_r.data());
    for (int ir = 0; ir < nrxx; ++ir)
    {
        dvxc_r[ir] -= dh_r[ir];
    }
}

} // namespace ModuleDFPT

#endif // __LIBXC
