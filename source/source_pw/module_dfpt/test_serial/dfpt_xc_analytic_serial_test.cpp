#include <gtest/gtest.h>

#ifdef __LIBXC

#include "source_pw/module_dfpt/dfpt_xc_analytic.h"

#include "source_base/constants.h"
#include "source_basis/module_pw/pw_basis.h"

#include <xc.h>

#include <cmath>
#include <complex>
#include <vector>

namespace
{

// simple serial PW_Basis on a cubic cell (rho-grid only); PW_Basis owns raw
// FFT state, so it is initialized in place instead of returned by value
void make_basis(ModulePW::PW_Basis& pw, const double lat0, const double a, const double ecut)
{
    const ModuleBase::Matrix3 latvec(a, 0.0, 0.0, 0.0, a, 0.0, 0.0, 0.0, a);
    pw.initgrids(lat0, latvec, ecut);
    pw.initparameters(false, ecut);
    pw.fft_bundle.initfftmode(0);
    pw.setuptransform();
    pw.collect_local_pw();
}

// smooth positive density: gaussian bump at the cell center
std::vector<double> gaussian_rho(const ModulePW::PW_Basis& pw, const double amp)
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
std::vector<double> cosine_amplitude(const ModulePW::PW_Basis& pw, const double amp)
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
void fft_gradient(const ModulePW::PW_Basis& pw, const std::vector<double>& f, std::vector<std::vector<double>>& grad)
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
void fft_divergence(const ModulePW::PW_Basis& pw,
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
void libxc_vxc(const std::vector<int>& func_id,
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

} // namespace

TEST(DFPTXcAnalytic, LdaKernelMatchesVrhofiniteDifference)
{
    ModulePW::PW_Basis pw;
    make_basis(pw, 1.8897261255, 6.0, 20.0);
    const std::vector<double> rho0 = gaussian_rho(pw, 0.3);
    const std::vector<int> func_id = {XC_LDA_X, XC_LDA_C_PZ};
    ModuleDFPT::XC_First_Order_Analytic kernel(func_id, 1, rho0, &pw);

    // constant amplitude A = 1: dvxc[ir] must equal k_rr[ir]
    std::vector<std::complex<double>> a(pw.nrxx, std::complex<double>(1.0, 0.0));
    std::vector<std::complex<double>> dvxc;
    kernel.apply(a, dvxc);

    // reference: central difference of the LDA potential. libxc vxc is
    // already d(rho e)/d rho, so the kernel is its plain derivative
    const double eta = 1.0e-5;
    std::vector<double> rho_p(pw.nrxx);
    std::vector<double> rho_m(pw.nrxx);
    for (int ir = 0; ir < pw.nrxx; ++ir)
    {
        rho_p[ir] = rho0[ir] + eta;
        rho_m[ir] = rho0[ir] - eta;
    }
    std::vector<double> v_p(pw.nrxx, 0.0);
    std::vector<double> v_m(pw.nrxx, 0.0);
    for (size_t ifunc = 0; ifunc < func_id.size(); ++ifunc)
    {
        xc_func_type func;
        ASSERT_EQ(xc_func_init(&func, func_id[ifunc], XC_UNPOLARIZED), 0);
        std::vector<double> zk(pw.nrxx, 0.0);
        std::vector<double> vxc(pw.nrxx, 0.0);
        xc_lda_exc_vxc(&func, pw.nrxx, rho_p.data(), zk.data(), vxc.data());
        for (int ir = 0; ir < pw.nrxx; ++ir)
        {
            v_p[ir] += vxc[ir];
        }
        xc_lda_exc_vxc(&func, pw.nrxx, rho_m.data(), zk.data(), vxc.data());
        for (int ir = 0; ir < pw.nrxx; ++ir)
        {
            v_m[ir] += vxc[ir];
        }
        xc_func_end(&func);
    }
    // sigma = |grad rho0|^2: the FD reference has no kernel guards, so the
    // comparison is restricted to the region where the analytic kernel is
    // active (rho > epsr, sigma > epsg) and dense enough that the FD
    // truncation error (grows like rho^-7/3 in the tail) stays controlled
    std::vector<std::vector<double>> grho0_l;
    fft_gradient(pw, rho0, grho0_l);
    double max_abs = 0.0;
    double max_err = 0.0;
    int n_checked = 0;
    for (int ir = 0; ir < pw.nrxx; ++ir)
    {
        if (rho0[ir] < 0.01)
        {
            continue;
        }
        double s2 = 0.0;
        for (int ipol = 0; ipol < 3; ++ipol)
        {
            s2 += grho0_l[ipol][ir] * grho0_l[ipol][ir];
        }
        if (s2 < 1.0e-8)
        {
            continue;
        }
        const double ref = 2.0 * (v_p[ir] - v_m[ir]) / (2.0 * eta); // Ry
        max_abs = std::max(max_abs, std::abs(ref));
        max_err = std::max(max_err, std::abs(ref - dvxc[ir].real()));
        ++n_checked;
    }
    EXPECT_GT(n_checked, pw.nrxx / 10); // the checked region is representative
    EXPECT_GT(max_abs, 1.0); // nontrivial kernel scale (Ry bohr^3)
    EXPECT_LT(max_err / max_abs, 1.0e-6);
}

TEST(DFPTXcAnalytic, GgaQ0MatchesFieldPerturbationFiniteDifference)
{
    ModulePW::PW_Basis pw;
    make_basis(pw, 1.8897261255, 6.0, 20.0);
    const std::vector<double> rho0 = gaussian_rho(pw, 0.3);
    const std::vector<double> amp = cosine_amplitude(pw, 0.02);
    const std::vector<int> func_id = {XC_GGA_X_PBE, XC_GGA_C_PBE};
    ModuleDFPT::XC_First_Order_Analytic kernel(func_id, 2, rho0, &pw);
    const ModuleBase::Vector3<double> q0(0.0, 0.0, 0.0);
    kernel.set_context(q0);

    std::vector<std::complex<double>> a(pw.nrxx);
    for (int ir = 0; ir < pw.nrxx; ++ir)
    {
        a[ir] = std::complex<double>(amp[ir], 0.0);
    }
    std::vector<std::complex<double>> dvxc;
    kernel.apply(a, dvxc);

    // reference: central difference of the full GGA potential
    //   v_xc = vrho - div(2 vsigma grad rho)
    // between rho0 +- eta A (Hartree -> Ry by the 2.0 below)
    std::vector<std::vector<double>> grho0;
    fft_gradient(pw, rho0, grho0);
    std::vector<std::vector<double>> gamp;
    fft_gradient(pw, amp, gamp);
    const double eta = 2.0e-6;
    std::vector<double> rho_p(pw.nrxx);
    std::vector<double> rho_m(pw.nrxx);
    std::vector<double> sig_p(pw.nrxx);
    std::vector<double> sig_m(pw.nrxx);
    for (int ir = 0; ir < pw.nrxx; ++ir)
    {
        rho_p[ir] = rho0[ir] + eta * amp[ir];
        rho_m[ir] = rho0[ir] - eta * amp[ir];
        double sp = 0.0;
        double sm = 0.0;
        for (int ipol = 0; ipol < 3; ++ipol)
        {
            const double gp = grho0[ipol][ir] + eta * gamp[ipol][ir];
            const double gm = grho0[ipol][ir] - eta * gamp[ipol][ir];
            sp += gp * gp;
            sm += gm * gm;
        }
        sig_p[ir] = sp;
        sig_m[ir] = sm;
    }
    std::vector<double> vr_p;
    std::vector<double> vs_p;
    std::vector<double> vr_m;
    std::vector<double> vs_m;
    libxc_vxc(func_id, rho_p, sig_p, vr_p, vs_p);
    libxc_vxc(func_id, rho_m, sig_m, vr_m, vs_m);
    std::vector<std::vector<double>> h_p(3, std::vector<double>(pw.nrxx, 0.0));
    std::vector<std::vector<double>> h_m(3, std::vector<double>(pw.nrxx, 0.0));
    for (int ir = 0; ir < pw.nrxx; ++ir)
    {
        for (int ipol = 0; ipol < 3; ++ipol)
        {
            const double gp = grho0[ipol][ir] + eta * gamp[ipol][ir];
            const double gm = grho0[ipol][ir] - eta * gamp[ipol][ir];
            h_p[ipol][ir] = 2.0 * vs_p[ir] * gp;
            h_m[ipol][ir] = 2.0 * vs_m[ir] * gm;
        }
    }
    std::vector<double> div_p;
    std::vector<double> div_m;
    fft_divergence(pw, h_p, div_p);
    fft_divergence(pw, h_m, div_m);
    // same region restriction as the LDA test: compare where the kernel is
    // active and the FD truncation is controlled
    double l2_num = 0.0;
    double l2_den = 0.0;
    double max_abs = 0.0;
    int n_checked = 0;
    for (int ir = 0; ir < pw.nrxx; ++ir)
    {
        if (rho0[ir] < 0.01)
        {
            continue;
        }
        double s2 = 0.0;
        for (int ipol = 0; ipol < 3; ++ipol)
        {
            s2 += grho0[ipol][ir] * grho0[ipol][ir];
        }
        if (s2 < 1.0e-8)
        {
            continue;
        }
        const double v_p = 2.0 * (vr_p[ir] - div_p[ir]);
        const double v_m = 2.0 * (vr_m[ir] - div_m[ir]);
        const double ref = (v_p - v_m) / (2.0 * eta);
        const double err = ref - dvxc[ir].real();
        l2_num += err * err;
        l2_den += ref * ref;
        max_abs = std::max(max_abs, std::abs(ref));
        ++n_checked;
    }
    EXPECT_GT(n_checked, pw.nrxx / 10); // the checked region is representative
    EXPECT_GT(max_abs, 1.0e-3); // nontrivial response scale
    // the eta^2 FD truncation of the reference bounds the agreement around
    // 1e-6..1e-5; structural kernel/assembly errors show up at percent level
    EXPECT_LT(std::sqrt(l2_num / l2_den), 1.0e-5);
}

TEST(DFPTXcAnalytic, CarrierTermAppearsAtFiniteQOnlyForGga)
{
    ModulePW::PW_Basis pw;
    make_basis(pw, 1.8897261255, 6.0, 20.0);
    const std::vector<double> rho0 = gaussian_rho(pw, 0.3);
    const std::vector<double> amp = cosine_amplitude(pw, 0.02);
    std::vector<std::complex<double>> a(pw.nrxx);
    for (int ir = 0; ir < pw.nrxx; ++ir)
    {
        a[ir] = std::complex<double>(amp[ir], 0.3 * amp[ir]);
    }

    // GGA: the finite-q response must differ from the Gamma response
    // (the i(q.f_grad)A carrier term); LDA must stay q-independent
    const std::vector<int> gga_id = {XC_GGA_X_PBE, XC_GGA_C_PBE};
    ModuleDFPT::XC_First_Order_Analytic gga(gga_id, 2, rho0, &pw);
    const ModuleBase::Vector3<double> q0(0.0, 0.0, 0.0);
    gga.set_context(q0);
    std::vector<std::complex<double>> d0;
    gga.apply(a, d0);
    const ModuleBase::Vector3<double> qx(0.25, 0.0, 0.0);
    gga.set_context(qx);
    std::vector<std::complex<double>> dq;
    gga.apply(a, dq);
    double diff = 0.0;
    double scale = 0.0;
    for (int ir = 0; ir < pw.nrxx; ++ir)
    {
        diff = std::max(diff, std::abs(dq[ir] - d0[ir]));
        scale = std::max(scale, std::abs(d0[ir]));
    }
    EXPECT_GT(diff / scale, 1.0e-3);

    const std::vector<int> lda_id = {XC_LDA_X, XC_LDA_C_PZ};
    ModuleDFPT::XC_First_Order_Analytic lda(lda_id, 1, rho0, &pw);
    lda.set_context(q0);
    std::vector<std::complex<double>> l0;
    lda.apply(a, l0);
    lda.set_context(qx);
    std::vector<std::complex<double>> lq;
    lda.apply(a, lq);
    double ldiff = 0.0;
    for (int ir = 0; ir < pw.nrxx; ++ir)
    {
        ldiff = std::max(ldiff, std::abs(lq[ir] - l0[ir]));
    }
    EXPECT_EQ(ldiff, 0.0);
}

#endif // __LIBXC
