#include "xc_functional_gga_noncol_sf_builtin.h"

#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_base/vector3.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/unitcell.h"
#include "source_io/module_parameter/parameter.h"
#include "xc_functional.h"

#include <cmath>
#include <vector>

namespace ModuleXC
{
namespace NCGGA_SF_Builtin
{

std::tuple<double, double, ModuleBase::matrix> v_xc_ncgga_sf_builtin(
    const int& nrxx, const double& omega, const double tpiba, const Charge* const chr)
{
    ModuleBase::TITLE("XC_Functional", "v_xc_ncgga_sf_builtin");
    ModuleBase::timer::start("XC_Functional", "v_xc_ncgga_sf_builtin");

    if (PARAM.inp.nspin != 4 || (!PARAM.globalv.domag && !PARAM.globalv.domag_z))
        throw std::domain_error("v_xc_ncgga_sf_builtin requires NSPIN==4.");

    // ======================================================================
    // Scalmani-Frisch (SF) builtin for noncollinear GGA (gga_grad=3)
    //
    // Reference: Scalmani & Frisch, JCTC 8, 1069 (2012)
    //
    // Key formulas:
    //   rho_up   = 0.5*(rho + |m|),    rho_dn = 0.5*(rho - |m|)
    //   m_hat    = m / |m|             (magnetization unit vector)
    //
    // Gradient decomposition (chain rule via m_hat):
    //   grad(rho_up) = 0.5 grad(rho) + 0.5 m_hat_mu * grad(m_mu)
    //   grad(rho_dn) = 0.5 grad(rho) - 0.5 m_hat_mu * grad(m_mu)
    //
    // For GGA, the sigma = |grad(rho_up)|^2 etc. are computed from gdr1,gdr2,
    // then passed to spin-polarized XC functionals (xc_spin, gcx_spin, gcc_spin).
    // The resulting potential is converted back to nspin=4 representation:
    //   v_tot = 0.5*(v_up + v_dn)
    //   v_mag = 0.5*(v_up - v_dn) * m_hat
    // ======================================================================

    ModulePW::PW_Basis* rhopw = chr->rhopw;
    const int npw = rhopw->npw;
    const double e2 = ModuleBase::e2;
    constexpr double vanishing = 1e-10;
    constexpr double epsr = 1e-6;
    const double fac = 0.5;
    const bool is_gga = (XC_Functional::get_func_type() == 2 || XC_Functional::get_func_type() == 4);

    // Step 1: decompose charge density into spin-up/spin-down and compute m_hat
    std::vector<double> rhotmp1(nrxx), rhotmp2(nrxx), amag(nrxx);
    std::vector<double> mag_part(3 * nrxx, 0.0);

    for (int ir = 0; ir < nrxx; ++ir)
    {
        const double mx = chr->rho[1][ir], my = chr->rho[2][ir], mz = chr->rho[3][ir];
        amag[ir] = std::sqrt(mx * mx + my * my + mz * mz);
        rhotmp1[ir] = 0.5 * (chr->rho[0][ir] + amag[ir]);
        rhotmp2[ir] = 0.5 * (chr->rho[0][ir] - amag[ir]);
        if (amag[ir] > 1e-12)
        {
            mag_part[ir] = mx / amag[ir];
            mag_part[ir + nrxx] = my / amag[ir];
            mag_part[ir + 2 * nrxx] = mz / amag[ir];
        }
    }
    for (int ir = 0; ir < nrxx; ++ir)
    {
        rhotmp1[ir] += fac * chr->rho_core[ir];
        rhotmp2[ir] += fac * chr->rho_core[ir];
    }

    std::vector<std::complex<double>> rhogsum1(npw), tmp_recip(npw);
    rhopw->real2recip(chr->rho[0], rhogsum1.data());
    for (int ig = 0; ig < npw; ++ig)
        rhogsum1[ig] += chr->rhog_core[ig];

    std::vector<ModuleBase::Vector3<double>> gdr1(nrxx), gdr2(nrxx);
    std::vector<ModuleBase::Vector3<double>> gdr_mag(nrxx);
    XC_Functional::grad_rho(rhogsum1.data(), gdr1.data(), rhopw, tpiba);

    for (int ir = 0; ir < nrxx; ++ir)
    {
        gdr_mag[ir] = gdr1[ir];
        gdr1[ir] = 0.5 * gdr_mag[ir];
        gdr2[ir] = 0.5 * gdr_mag[ir];
    }
    for (int is = 1; is <= 3; ++is)
    {
        rhopw->real2recip(chr->rho[is], tmp_recip.data());
        XC_Functional::grad_rho(tmp_recip.data(), gdr_mag.data(), rhopw, tpiba);
        const double* mp = mag_part.data() + (is - 1) * nrxx;
        for (int ir = 0; ir < nrxx; ++ir)
        {
            const ModuleBase::Vector3<double> g = 0.5 * gdr_mag[ir] * mp[ir];
            gdr1[ir] += g;
            gdr2[ir] -= g;
        }
    }

    // Step 3: LDA contribution (xc_spin) and convert to nspin=4 potential
    //   v(0,:) = 0.5*(vxc[0] + vxc[1])    -- total density channel
    //   v(1..3,:) = 0.5*(vxc[0] - vxc[1]) * m_mu / |m|    -- magnetization channels
    double etxc = 0, vtxc = 0;
    ModuleBase::matrix v(PARAM.inp.nspin, nrxx);

    for (int ir = 0; ir < nrxx; ++ir)
    {
        const double arho = std::abs(chr->rho[0][ir] + chr->rho_core[ir]);
        if (arho <= vanishing)
            continue;

        double zeta = amag[ir] / arho;
        if (std::abs(zeta) > 1.0)
            zeta = (zeta > 0) ? 1.0 : -1.0;
        double exc = 0, vxc[2] = {0, 0};
        XC_Functional::xc_spin(arho, zeta, exc, vxc[0], vxc[1]);

        v(0, ir) = e2 * 0.5 * (vxc[0] + vxc[1]);
        vtxc += v(0, ir) * chr->rho[0][ir];

        if (amag[ir] > vanishing)
        {
            const double vs = e2 * 0.5 * (vxc[0] - vxc[1]);
            const double inv_a = 1.0 / amag[ir];
            for (int mu = 1; mu < 4; ++mu)
            {
                v(mu, ir) = vs * chr->rho[mu][ir] * inv_a;
                vtxc += v(mu, ir) * chr->rho[mu][ir];
            }
        }
        etxc += e2 * exc * arho;
    }

    // Step 4: GGA contribution (gcx_spin + gcc_spin) and divergence correction (grad_dot)
    //   The gradient correction uses the same SF-decomposed gdr1, gdr2 and converts
    //   back to nspin=4 via m_hat, matching the LDA conversion above.
    if (is_gga)
    {
        double etxcgc = 0, vtxcgc = 0;
        std::vector<double> vup_gga(nrxx, 0), vdw_gga(nrxx, 0);
        std::vector<ModuleBase::Vector3<double>> h1(nrxx), h2(nrxx);

        for (int ir = 0; ir < nrxx; ++ir)
        {
            double sx = 0, v1xup = 0, v1xdw = 0, v2xup = 0, v2xdw = 0;
            double sc = 0, v1cup = 0, v1cdw = 0, v2c = 0;
            const double grho2a = gdr1[ir] * gdr1[ir];
            const double grho2b = gdr2[ir] * gdr2[ir];
            const double rh = rhotmp1[ir] + rhotmp2[ir];

            XC_Functional::gcx_spin(rhotmp1[ir], rhotmp2[ir], grho2a, grho2b,
                                    sx, v1xup, v1xdw, v2xup, v2xdw);

            if (rh > epsr)
            {
                double zeta = (rhotmp1[ir] - rhotmp2[ir]) / rh;
                zeta = std::fabs(zeta);
                const double grh2 = (gdr1[ir] + gdr2[ir]) * (gdr1[ir] + gdr2[ir]);
                XC_Functional::gcc_spin(rh, zeta, grh2, sc, v1cup, v1cdw, v2c);
            }

            vup_gga[ir] = e2 * (v1xup + v1cup);
            vdw_gga[ir] = e2 * (v1xdw + v1cdw);

            const double v2cup = v2c, v2cdw = v2c, v2cud = v2c;
            h1[ir] = e2 * ((v2xup + v2cup) * gdr1[ir] + v2cud * gdr2[ir]);
            h2[ir] = e2 * ((v2xdw + v2cdw) * gdr2[ir] + v2cud * gdr1[ir]);

            vtxcgc += vup_gga[ir] * (rhotmp1[ir] - chr->rho_core[ir] * fac);
            vtxcgc += vdw_gga[ir] * (rhotmp2[ir] - chr->rho_core[ir] * fac);
            etxcgc += e2 * (sx + sc);
        }

        for (int ir = 0; ir < nrxx; ++ir)
        {
            v(0, ir) += 0.5 * (vup_gga[ir] + vdw_gga[ir]);
            const double vdiff = 0.5 * (vup_gga[ir] - vdw_gga[ir]);
            for (int mu = 1; mu < 4; ++mu)
            {
                v(mu, ir) += vdiff * mag_part[ir + (mu - 1) * nrxx];
            }
        }

        std::vector<double> dh(nrxx);
        std::vector<ModuleBase::Vector3<double>> tmp_h(nrxx);

        for (int ir = 0; ir < nrxx; ++ir)
            tmp_h[ir] = 0.5 * (h1[ir] + h2[ir]);
        XC_Functional::grad_dot(tmp_h.data(), dh.data(), rhopw, tpiba);
        for (int ir = 0; ir < nrxx; ++ir)
            v(0, ir) -= dh[ir];
        double sum = 0;
        for (int ir = 0; ir < nrxx; ++ir)
            sum += dh[ir] * chr->rho[0][ir];
        vtxcgc -= sum;

        for (int mu = 1; mu < 4; ++mu)
        {
            const double* mp = mag_part.data() + (mu - 1) * nrxx;
            for (int ir = 0; ir < nrxx; ++ir)
                tmp_h[ir] = 0.5 * (h1[ir] - h2[ir]) * mp[ir];
            XC_Functional::grad_dot(tmp_h.data(), dh.data(), rhopw, tpiba);
            for (int ir = 0; ir < nrxx; ++ir)
                v(mu, ir) -= dh[ir];
            double sum_mu = 0;
            for (int ir = 0; ir < nrxx; ++ir)
                sum_mu += dh[ir] * chr->rho[mu][ir];
            vtxcgc -= sum_mu;
        }

        etxc += etxcgc;
        vtxc += vtxcgc;
    }

#ifdef __MPI
    Parallel_Reduce::reduce_pool(etxc);
    Parallel_Reduce::reduce_pool(vtxc);
#endif
    etxc *= omega / rhopw->nxyz;
    vtxc *= omega / rhopw->nxyz;

    ModuleBase::timer::end("XC_Functional", "v_xc_ncgga_sf_builtin");
    return std::make_tuple(etxc, vtxc, std::move(v));
}

void gradcorr_ncgga_sf_builtin(const Charge* const chr, ModulePW::PW_Basis* rhopw,
                                const UnitCell* ucell, std::vector<double>& stress_gga)
{
    stress_gga.assign(9, 0.0);

    const int nrxx = rhopw->nrxx;
    const int npw = rhopw->npw;
    const double e2 = ModuleBase::e2;
    constexpr double epsr = 1.0e-6;
    constexpr double small = 1.0e-10;
    const double fac = 0.5;

    std::vector<double> rhotmp1(nrxx), rhotmp2(nrxx), amag_arr(nrxx);
    std::vector<double> mag_part(3 * nrxx, 0.0);

    for (int ir = 0; ir < nrxx; ++ir)
    {
        const double mx = chr->rho[1][ir], my = chr->rho[2][ir], mz = chr->rho[3][ir];
        double amag = std::sqrt(mx * mx + my * my + mz * mz);
        amag_arr[ir] = amag;
        rhotmp1[ir] = 0.5 * (chr->rho[0][ir] + amag);
        rhotmp2[ir] = 0.5 * (chr->rho[0][ir] - amag);
        if (amag > 1e-12)
        {
            mag_part[ir] = mx / amag;
            mag_part[ir + nrxx] = my / amag;
            mag_part[ir + 2 * nrxx] = mz / amag;
        }
    }
    for (int ir = 0; ir < nrxx; ++ir)
    {
        rhotmp1[ir] += fac * chr->rho_core[ir];
        rhotmp2[ir] += fac * chr->rho_core[ir];
    }

    std::vector<std::complex<double>> rhogsum1(npw), tmp_recip(npw);
    rhopw->real2recip(chr->rho[0], rhogsum1.data());
    for (int ig = 0; ig < npw; ++ig)
        rhogsum1[ig] += chr->rhog_core[ig];

    std::vector<ModuleBase::Vector3<double>> gdr1(nrxx), gdr2(nrxx);
    std::vector<ModuleBase::Vector3<double>> gdr_mag(nrxx);
    XC_Functional::grad_rho(rhogsum1.data(), gdr1.data(), rhopw, ucell->tpiba);

    for (int ir = 0; ir < nrxx; ++ir)
    {
        gdr_mag[ir] = gdr1[ir];
        gdr1[ir] = 0.5 * gdr_mag[ir];
        gdr2[ir] = 0.5 * gdr_mag[ir];
    }
    for (int is = 1; is <= 3; ++is)
    {
        rhopw->real2recip(chr->rho[is], tmp_recip.data());
        XC_Functional::grad_rho(tmp_recip.data(), gdr_mag.data(), rhopw, ucell->tpiba);
        const double* mp = mag_part.data() + (is - 1) * nrxx;
        for (int ir = 0; ir < nrxx; ++ir)
        {
            const ModuleBase::Vector3<double> g = 0.5 * gdr_mag[ir] * mp[ir];
            gdr1[ir] += g;
            gdr2[ir] -= g;
        }
    }

    for (int ir = 0; ir < nrxx; ++ir)
    {
        double sx = 0, v1xup = 0, v1xdw = 0, v2xup = 0, v2xdw = 0;
        double sc = 0, v1cup = 0, v1cdw = 0, v2c = 0;
        double v2cup = 0, v2cdw = 0, v2cud = 0;
        const double grho2a = gdr1[ir] * gdr1[ir];
        const double grho2b = gdr2[ir] * gdr2[ir];
        const double rh = rhotmp1[ir] + rhotmp2[ir];

        XC_Functional::gcx_spin(rhotmp1[ir], rhotmp2[ir], grho2a, grho2b,
                                sx, v1xup, v1xdw, v2xup, v2xdw);

        if (rh > epsr)
        {
            double zeta = (rhotmp1[ir] - rhotmp2[ir]) / rh;
            zeta = std::fabs(zeta);
            const double grh2 = (gdr1[ir] + gdr2[ir]) * (gdr1[ir] + gdr2[ir]);
            XC_Functional::gcc_spin(rh, zeta, grh2, sc, v1cup, v1cdw, v2c);
            v2cup = v2c;
            v2cdw = v2c;
            v2cud = v2c;
        }

        double tt1[3] = {gdr1[ir].x, gdr1[ir].y, gdr1[ir].z};
        double tt2[3] = {gdr2[ir].x, gdr2[ir].y, gdr2[ir].z};
        for (int l = 0; l < 3; l++)
        {
            for (int m = 0; m < l + 1; m++)
            {
                int ind = l * 3 + m;
                stress_gga[ind] += e2 * (tt1[l] * tt1[m] * v2xup + tt2[l] * tt2[m] * v2xdw);
                stress_gga[ind] += e2 * (tt1[l] * tt1[m] * v2cup + tt2[l] * tt2[m] * v2cdw
                                         + (tt1[l] * tt2[m] + tt2[l] * tt1[m]) * v2cud);
            }
        }
    }
}

} // namespace NCGGA_SF_Builtin
} // namespace ModuleXC
