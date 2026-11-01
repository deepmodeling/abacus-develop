#ifndef DFPT_XC_ANALYTIC_SERIAL_FIXTURE_H
#define DFPT_XC_ANALYTIC_SERIAL_FIXTURE_H

// Shared helper layer for the DFPT_XcAnalytic serial test (guarded by
// __LIBXC together with the test body): a cubic rho-grid basis builder,
// synthetic density/amplitude fields, FFT gradient/divergence references
// in physical units, and the summed libxc vrho/vsigma evaluation.

#ifdef __LIBXC

#include "source_base/constants.h"
#include "source_basis/module_pw/pw_basis.h"

#include "gtest/gtest.h"
#include <xc.h>
#include <xc_funcs.h>

#include <cmath>
#include <complex>
#include <vector>

namespace dfpt_xc_test_detail
{

// simple serial PW_Basis on a cubic cell (rho-grid only); PW_Basis owns raw
// FFT state, so it is initialized in place instead of returned by value
inline void make_basis(ModulePW::PW_Basis& pw, const double lat0, const double a, const double ecut)
{
    const ModuleBase::Matrix3 latvec(a, 0.0, 0.0, 0.0, a, 0.0, 0.0, 0.0, a);
    pw.initgrids(lat0, latvec, ecut);
    pw.initparameters(false, ecut);
    pw.fft_bundle.initfftmode(0);
    pw.setuptransform();
    pw.collect_local_pw();
}

// smooth positive density: gaussian bump at the cell center
inline std::vector<double> gaussian_rho(const ModulePW::PW_Basis& pw, const double amp)
{
    std::vector<double> rho(pw.nrxx, 0.0);
    for (int ix = 0; ix < pw.nx; ++ix)
    {
        for (int iy = 0; iy < pw.ny; ++iy)
        {
            for (int iz = 0; iz < pw.nz; ++iz)
            {
                const int ir = (ix * pw.ny + iy) * pw.nz + iz;
                const double x = static_cast<double>(ix) / pw.nx - 0.5;
                const double y = static_cast<double>(iy) / pw.ny - 0.5;
                const double z = static_cast<double>(iz) / pw.nz - 0.5;
                rho[ir] = amp * std::exp(-(x * x + y * y + z * z) / 0.08);
            }
        }
    }
    return rho;
}

// smooth real amplitude with gradient content
inline std::vector<double> cosine_amplitude(const ModulePW::PW_Basis& pw, const double amp)
{
    std::vector<double> a(pw.nrxx, 0.0);
    for (int ix = 0; ix < pw.nx; ++ix)
    {
        for (int iy = 0; iy < pw.ny; ++iy)
        {
            for (int iz = 0; iz < pw.nz; ++iz)
            {
                const int ir = (ix * pw.ny + iy) * pw.nz + iz;
                const double x = ModuleBase::TWO_PI * ix / pw.nx;
                const double y = ModuleBase::TWO_PI * iy / pw.ny;
                a[ir] = amp * (std::cos(x) + 0.5 * std::cos(y) + 0.3 * std::cos(x + y) * std::cos(2.0 * x));
            }
        }
    }
    return a;
}

// FFT gradient of a real field on the rho grid (physical units: pw gcar is
// dimensionless, the physical wavevector is tpiba * gcar)
inline void fft_gradient(const ModulePW::PW_Basis& pw,
                         const std::vector<double>& f,
                         std::vector<std::vector<double>>& grad)
{
    const int npw = pw.npw;
    std::vector<std::complex<double>> f_g(npw);
    pw.real2recip(f.data(), f_g.data());
    grad.assign(3, std::vector<double>(pw.nrxx, 0.0));
    std::vector<std::complex<double>> g_g(npw);
    for (int ipol = 0; ipol < 3; ++ipol)
    {
        for (int ig = 0; ig < npw; ++ig)
        {
            g_g[ig] = std::complex<double>(0.0, pw.tpiba * pw.gcar[ig][ipol]) * f_g[ig];
        }
        std::vector<std::complex<double>> g_r(pw.nrxx);
        pw.recip2real(g_g.data(), g_r.data());
        for (int ir = 0; ir < pw.nrxx; ++ir)
        {
            grad[ipol][ir] = g_r[ir].real();
        }
    }
}

// divergence of a real vector field on the rho grid (physical units)
inline void fft_divergence(const ModulePW::PW_Basis& pw,
                           const std::vector<std::vector<double>>& h,
                           std::vector<double>& div)
{
    const int npw = pw.npw;
    std::vector<std::complex<double>> dh_g(npw, std::complex<double>(0.0, 0.0));
    std::vector<std::complex<double>> h_g(npw);
    for (int ipol = 0; ipol < 3; ++ipol)
    {
        pw.real2recip(h[ipol].data(), h_g.data());
        for (int ig = 0; ig < npw; ++ig)
        {
            dh_g[ig] += std::complex<double>(0.0, pw.tpiba * pw.gcar[ig][ipol]) * h_g[ig];
        }
    }
    div.resize(pw.nrxx);
    std::vector<std::complex<double>> dh_r(pw.nrxx);
    pw.recip2real(dh_g.data(), dh_r.data());
    for (int ir = 0; ir < pw.nrxx; ++ir)
    {
        div[ir] = dh_r[ir].real();
    }
}

// sum of libxc vrho / vsigma over functional components (Hartree)
inline void libxc_vxc(const std::vector<int>& func_id,
                      const std::vector<double>& rho,
                      const std::vector<double>& sigma,
                      std::vector<double>& vrho,
                      std::vector<double>& vsigma)
{
    const int n = rho.size();
    vrho.assign(n, 0.0);
    vsigma.assign(n, 0.0);
    for (size_t ifunc = 0; ifunc < func_id.size(); ++ifunc)
    {
        xc_func_type func;
        ASSERT_EQ(xc_func_init(&func, func_id[ifunc], XC_UNPOLARIZED), 0);
        std::vector<double> vr(n, 0.0);
        std::vector<double> vs(n, 0.0);
        std::vector<double> zk(n, 0.0);
        xc_gga_exc_vxc(&func, n, rho.data(), sigma.data(), zk.data(), vr.data(), vs.data());
        for (int ir = 0; ir < n; ++ir)
        {
            vrho[ir] += vr[ir];
            vsigma[ir] += vs[ir];
        }
        xc_func_end(&func);
    }
}

} // namespace dfpt_xc_test_detail

#endif // __LIBXC

#endif // DFPT_XC_ANALYTIC_SERIAL_FIXTURE_H
