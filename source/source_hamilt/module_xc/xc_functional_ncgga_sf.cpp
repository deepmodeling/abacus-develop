#include "xc_functional_ncgga_sf.h"

#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_base/vector3.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_estate/module_charge/charge.h"
#include "xc_functional.h"
#include "xc_ncgga_radial.h"

#include <array>
#include <cmath>
#include <vector>

namespace ModuleXC
{
namespace NCGGA_SF_Builtin
{

std::tuple<double, double, ModuleBase::matrix> v_xc_ncgga_sf_builtin(const int& nrxx,
                                                                     const double& omega,
                                                                     const double tpiba,
                                                                     const Charge* const chr)
{
    ModuleBase::TITLE("XC_Functional", "v_xc_ncgga_sf_builtin");
    ModuleBase::timer::start("XC_Functional", "v_xc_ncgga_sf_builtin");

    // Regularized projected local-collinear energy: rho_s = N_s(n + rho_core, m),
    // g_s = sum_A J_sA G_h x_A. Reverse the same discrete graph, including
    // both -D_h(sum_s J_sA h_s) and the local Hessian response of N_s.
    ModulePW::PW_Basis* rhopw = chr->rhopw;
    const int npw = rhopw->npw;
    const double e2 = ModuleBase::e2;
    constexpr double vanishing = 1e-10;
    constexpr double epsr = 1e-6;
    const bool is_gga = (XC_Functional::get_func_type() == 2 || XC_Functional::get_func_type() == 4);

    std::vector<double> rhotmp1(nrxx);
    std::vector<double> rhotmp2(nrxx);
    std::vector<NcggaSpinMapPoint> spin_map(nrxx);

    for (int ir = 0; ir < nrxx; ++ir)
    {
        const double mx = chr->rho[1][ir];
        const double my = chr->rho[2][ir];
        const double mz = chr->rho[3][ir];

        const std::array<double, 3> magnetization = {{mx, my, mz}};
        const NcggaRadialPoint radial = make_ncgga_radial_point(magnetization, ncgga_lca_radial_eta());
        spin_map[ir] = make_ncgga_spin_map_point(chr->rho[0][ir] + chr->rho_core[ir], radial);
        rhotmp1[ir] = spin_map[ir].spin_density[0];
        rhotmp2[ir] = spin_map[ir].spin_density[1];
    }

    std::vector<std::complex<double>> rhogsum1(npw);
    std::vector<std::complex<double>> tmp_recip(npw);
    rhopw->real2recip(chr->rho[0], rhogsum1.data());
    for (int ig = 0; ig < npw; ++ig)
        rhogsum1[ig] += chr->rhog_core[ig];

    std::vector<ModuleBase::Vector3<double>> gdr1(nrxx);
    std::vector<ModuleBase::Vector3<double>> gdr2(nrxx);
    std::vector<ModuleBase::Vector3<double>> grad_rho(nrxx);
    std::array<std::vector<ModuleBase::Vector3<double>>, 3> grad_m;
    for (int mu = 0; mu < 3; ++mu)
    {
        grad_m[mu].resize(nrxx);
    }
    std::vector<ModuleBase::Vector3<double>> gdr_mag(nrxx);
    XC_Functional::grad_rho(rhogsum1.data(), gdr1.data(), rhopw, tpiba);

    for (int ir = 0; ir < nrxx; ++ir)
    {
        grad_rho[ir] = gdr1[ir];

        gdr1[ir] = spin_map[ir].jacobian(0, 0) * grad_rho[ir];
        gdr2[ir] = spin_map[ir].jacobian(1, 0) * grad_rho[ir];
    }
    for (int is = 1; is <= 3; ++is)
    {
        rhopw->real2recip(chr->rho[is], tmp_recip.data());
        XC_Functional::grad_rho(tmp_recip.data(), gdr_mag.data(), rhopw, tpiba);
        grad_m[is - 1] = gdr_mag;
        for (int ir = 0; ir < nrxx; ++ir)
        {

            gdr1[ir] += spin_map[ir].jacobian(0, is) * gdr_mag[ir];
            gdr2[ir] += spin_map[ir].jacobian(1, is) * gdr_mag[ir];
        }
    }

    double etxc = 0;
    double vtxc = 0;
    ModuleBase::matrix v(4, nrxx);

    for (int ir = 0; ir < nrxx; ++ir)
    {
        const double arho = spin_map[ir].absolute_density;
        if (arho <= vanishing)
            continue;

        double zeta = spin_map[ir].clipped_magnitude / arho;
        if (std::abs(zeta) > 1.0)
            zeta = (zeta > 0) ? 1.0 : -1.0;
        double exc = 0;
        double vxc[2] = {0, 0};
        XC_Functional::xc_spin(arho, zeta, exc, vxc[0], vxc[1]);

        for (int channel = 0; channel < 4; ++channel)
        {
            v(channel, ir)
                = e2 * (spin_map[ir].jacobian(0, channel) * vxc[0] + spin_map[ir].jacobian(1, channel) * vxc[1]);
        }

        etxc += e2 * exc * arho;
    }

    // Step 4: GGA contribution and variational divergence correction.
    if (is_gga)
    {
        double etxcgc = 0;
        std::vector<double> vup_gga(nrxx, 0);
        std::vector<double> vdw_gga(nrxx, 0);
        std::vector<ModuleBase::Vector3<double>> h1(nrxx);
        std::vector<ModuleBase::Vector3<double>> h2(nrxx);
        for (int ir = 0; ir < nrxx; ++ir)
        {
            double sx = 0;
            double v1xup = 0;
            double v1xdw = 0;
            double v2xup = 0;
            double v2xdw = 0;
            double sc = 0;
            double v1cup = 0;
            double v1cdw = 0;
            double v2c = 0;
            double grho2a = gdr1[ir] * gdr1[ir];
            double grho2b = gdr2[ir] * gdr2[ir];

            const double rh = rhotmp1[ir] + rhotmp2[ir];

            XC_Functional::gcx_spin(rhotmp1[ir], rhotmp2[ir], grho2a, grho2b, sx, v1xup, v1xdw, v2xup, v2xdw);

            if (rh > epsr)
            {
                const double zeta_input = std::fabs((rhotmp1[ir] - rhotmp2[ir]) / rh);
                double zeta = zeta_input;
                const double grh2 = (gdr1[ir] + gdr2[ir]) * (gdr1[ir] + gdr2[ir]);
                XC_Functional::gcc_spin(rh, zeta, grh2, sc, v1cup, v1cdw, v2c);
                if (zeta_input > 1.0 - epsr)
                {
                    // gcc_spin evaluates this branch at a fixed clipped zeta.
                    // Reverse that actual branch instead of differentiating
                    // through the discarded input polarization.
                    const double fixed_zeta_density_derivative = 0.5 * ((1.0 + zeta) * v1cup + (1.0 - zeta) * v1cdw);
                    v1cup = fixed_zeta_density_derivative;
                    v1cdw = fixed_zeta_density_derivative;
                }
            }

            vup_gga[ir] = e2 * (v1xup + v1cup);
            vdw_gga[ir] = e2 * (v1xdw + v1cdw);

            const double v2cup = v2c;
            const double v2cdw = v2c;
            const double v2cud = v2c;
            h1[ir] = e2 * ((v2xup + v2cup) * gdr1[ir] + v2cud * gdr2[ir]);
            h2[ir] = e2 * ((v2xdw + v2cdw) * gdr2[ir] + v2cud * gdr1[ir]);

            etxcgc += e2 * (sx + sc);
        }

        for (int ir = 0; ir < nrxx; ++ir)
        {

            for (int channel = 0; channel < 4; ++channel)
            {
                v(channel, ir) += spin_map[ir].jacobian(0, channel) * vup_gga[ir]
                                  + spin_map[ir].jacobian(1, channel) * vdw_gga[ir];
            }
        }

        std::vector<double> dh(nrxx);
        std::vector<ModuleBase::Vector3<double>> tmp_h(nrxx);

        // Exact reverse of g_s=sum_A J_sA G_h(x_A):
        //   v_B = -D_h(sum_s J_sB h_s)
        //         + sum_s,A dJ_sA/dx_B h_s.G_h(x_A).
        for (int channel = 0; channel < 4; ++channel)
        {
            for (int ir = 0; ir < nrxx; ++ir)
            {
                tmp_h[ir] = spin_map[ir].jacobian(0, channel) * h1[ir] + spin_map[ir].jacobian(1, channel) * h2[ir];
            }
            XC_Functional::grad_dot(tmp_h.data(), dh.data(), rhopw, tpiba);
            for (int ir = 0; ir < nrxx; ++ir)
            {
                v(channel, ir) -= dh[ir];
                if (channel == 0 || spin_map[ir].saturated)
                {
                    continue;
                }
                const ModuleBase::Vector3<double> spin_flux = 0.5 * (h1[ir] - h2[ir]);
                double local_response = 0.0;
                for (int nu = 0; nu < 3; ++nu)
                {
                    local_response += spin_map[ir].radial.jacobian(nu, channel - 1) * (spin_flux * grad_m[nu][ir]);
                }
                v(channel, ir) += local_response;
            }
        }

        etxc += etxcgc;
    }

    // vtxc uses the same completed four-component potential returned to the
    // caller.  This unifies the bookkeeping for both modes.
    vtxc = 0.0;
    for (int ir = 0; ir < nrxx; ++ir)
    {
        for (int is = 0; is < 4; ++is)
        {
            vtxc += v(is, ir) * chr->rho[is][ir];
        }
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

} // namespace NCGGA_SF_Builtin
} // namespace ModuleXC
