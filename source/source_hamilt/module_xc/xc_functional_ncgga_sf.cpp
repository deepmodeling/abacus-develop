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

    // Mapped local-collinear energy: rho_s = N_s(n + rho_core, m),
    // g_s = G_h(rho_s). Reverse the FFT gradient before the local Jacobian:
    // v_A = sum_s J_sA [f_rho_s - D_h(h_s)]. No discrete product rule is used.
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

    std::vector<std::complex<double>> reciprocal(npw);
    std::vector<ModuleBase::Vector3<double>> gdr1(nrxx);
    std::vector<ModuleBase::Vector3<double>> gdr2(nrxx);
    rhopw->real2recip(rhotmp1.data(), reciprocal.data());
    XC_Functional::grad_rho(reciprocal.data(), gdr1.data(), rhopw, tpiba);
    rhopw->real2recip(rhotmp2.data(), reciprocal.data());
    XC_Functional::grad_rho(reciprocal.data(), gdr2.data(), rhopw, tpiba);

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

        std::vector<double> divergence_up(nrxx);
        std::vector<double> divergence_down(nrxx);
        XC_Functional::grad_dot(h1.data(), divergence_up.data(), rhopw, tpiba);
        XC_Functional::grad_dot(h2.data(), divergence_down.data(), rhopw, tpiba);
        for (int ir = 0; ir < nrxx; ++ir)
        {
            for (int channel = 0; channel < 4; ++channel)
            {
                v(channel, ir) -= spin_map[ir].jacobian(0, channel) * divergence_up[ir]
                                  + spin_map[ir].jacobian(1, channel) * divergence_down[ir];
            }
        }

        etxc += etxcgc;
    }

    // vtxc uses the same completed four-component potential returned to the
    // caller, including all four density channels.
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

void gradcorr_ncgga_lca_builtin(const Charge* const chr,
                                ModulePW::PW_Basis* rhopw,
                                const double tpiba,
                                std::vector<double>& stress_gga)
{
    stress_gga.assign(9, 0.0);

    const int nrxx = rhopw->nrxx;
    const int npw = rhopw->npw;
    const double e2 = ModuleBase::e2;
    constexpr double epsr = 1.0e-6;

    // Use exactly the spin-density gradients in the discrete XC energy.
    // For a metric variation, fixed real-grid spin densities remain fixed.
    std::vector<NcggaSpinMapPoint> spin_map(nrxx);
    for (int ir = 0; ir < nrxx; ++ir)
    {
        const std::array<double, 3> magnetization
            = {{chr->rho[1][ir], chr->rho[2][ir], chr->rho[3][ir]}};
        const double eta = ncgga_lca_radial_eta();
        const NcggaRadialPoint radial = make_ncgga_radial_point(magnetization, eta);
        const double total_density = chr->rho[0][ir] + chr->rho_core[ir];
        spin_map[ir] = make_ncgga_spin_map_point(total_density, radial);
    }
    std::array<std::vector<ModuleBase::Vector3<double>>, 2> spin_gradient;
    std::vector<double> spin_density(nrxx);
    std::vector<std::complex<double>> reciprocal(npw);
    for (int spin = 0; spin < 2; ++spin)
    {
        spin_gradient[spin].resize(nrxx);
        for (int ir = 0; ir < nrxx; ++ir)
        {
            spin_density[ir] = spin_map[ir].spin_density[spin];
        }
        rhopw->real2recip(spin_density.data(), reciprocal.data());
        XC_Functional::grad_rho(reciprocal.data(), spin_gradient[spin].data(), rhopw, tpiba);
    }

    for (int ir = 0; ir < nrxx; ++ir)
    {
        const double rho_up = spin_map[ir].spin_density[0];
        const double rho_down = spin_map[ir].spin_density[1];
        const ModuleBase::Vector3<double>& grad_up = spin_gradient[0][ir];
        const ModuleBase::Vector3<double>& grad_down = spin_gradient[1][ir];

        double sx = 0.0;
        double v1xup = 0.0;
        double v1xdw = 0.0;
        double v2xup = 0.0;
        double v2xdw = 0.0;
        XC_Functional::gcx_spin(rho_up,
                                rho_down,
                                grad_up * grad_up,
                                grad_down * grad_down,
                                sx,
                                v1xup,
                                v1xdw,
                                v2xup,
                                v2xdw);

        double sc = 0.0;
        double v1cup = 0.0;
        double v1cdw = 0.0;
        double v2c = 0.0;
        const double rho = rho_up + rho_down;
        if (rho > epsr)
        {
            double zeta = std::fabs((rho_up - rho_down) / rho);
            const ModuleBase::Vector3<double> grad_rho = grad_up + grad_down;
            XC_Functional::gcc_spin(rho, zeta, grad_rho * grad_rho, sc, v1cup, v1cdw, v2c);
        }

        const ModuleBase::Vector3<double> h_up = e2 * ((v2xup + v2c) * grad_up + v2c * grad_down);
        const ModuleBase::Vector3<double> h_down = e2 * ((v2xdw + v2c) * grad_down + v2c * grad_up);
        const double grad_up_component[3] = {grad_up.x, grad_up.y, grad_up.z};
        const double grad_down_component[3] = {grad_down.x, grad_down.y, grad_down.z};
        const double h_up_component[3] = {h_up.x, h_up.y, h_up.z};
        const double h_down_component[3] = {h_down.x, h_down.y, h_down.z};
        for (int row = 0; row < 3; ++row)
        {
            for (int column = 0; column <= row; ++column)
            {
                stress_gga[row * 3 + column] += h_up_component[row] * grad_up_component[column]
                                                + h_down_component[row] * grad_down_component[column];
            }
        }
    }
}

} // namespace NCGGA_SF_Builtin
} // namespace ModuleXC
