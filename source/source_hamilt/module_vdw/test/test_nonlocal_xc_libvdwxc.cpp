#include "../nonlocal_xc_libvdwxc.h"

#include "source_basis/module_pw/pw_basis.h"
#include "source_hamilt/module_xc/xc_functional.h"

#include <complex.h>
#include <vdwxc.h>

#include <cmath>
#include <complex>
#include <gtest/gtest.h>
#include <vector>

namespace
{
const double pi = std::acos(-1.0);

class LibvdwxcAdapterTest : public testing::Test
{
  protected:
    ModulePW::PW_Basis pw{"cpu", "double"};

    void SetUp() override
    {
        // The oblique cell catches row/column-major cell mapping mistakes.
        const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0, -0.5, std::sqrt(0.75), 0.0, 0.0, 0.0, 2.0);
        pw.initgrids(10.0, lattice, 6, 8, 12);
        pw.initparameters(false, 4.0, 1, true);
        pw.setuptransform();
        pw.collect_local_pw();
    }

    std::vector<double> density(double offset, double amplitude) const
    {
        std::vector<double> values(pw.nrxx);
        for (int x = 0; x < pw.nx; ++x)
            for (int y = 0; y < pw.ny; ++y)
                for (int z = 0; z < pw.nz; ++z)
                {
                    const int index = (x * pw.ny + y) * pw.nz + z;
                    const double ax = 2.0 * pi * x / pw.nx;
                    const double ay = 2.0 * pi * y / pw.ny;
                    const double az = 2.0 * pi * z / pw.nz;
                    values[index] = offset + amplitude * (1.0 + 0.25 * std::sin(ax) + 0.15 * std::cos(ay + az));
                }
        return values;
    }
};

#ifdef __LIBVDWXC
TEST_F(LibvdwxcAdapterTest, MatchesDirectSerialApiForNonuniformDensity)
{
    const double b = 6.3;
    const double c = 0.0093;
    const std::vector<double> total = density(0.004, 0.006);
    const std::vector<double> valence = density(0.002, 0.003);
    vdw::LibvdwxcNonlocalXC adapter(pw, b, c);
    const vdw::NonlocalXCResult adapted = adapter.evaluate(pw, total, valence);

    vdwxc_data direct = vdwxc_new(VDWXC_RVV10);
    ASSERT_NE(direct, nullptr);
    ASSERT_EQ(vdwxc_set_rvv10_parameters(direct, b, c), 0);
    const auto& a = pw.latvec;
    vdwxc_set_unit_cell(direct, pw.nx, pw.ny, pw.nz,
                        a.e11 * pw.lat0, a.e21 * pw.lat0, a.e31 * pw.lat0,
                        a.e12 * pw.lat0, a.e22 * pw.lat0, a.e32 * pw.lat0,
                        a.e13 * pw.lat0, a.e23 * pw.lat0, a.e33 * pw.lat0);
    vdwxc_init_serial(direct);

    std::vector<std::complex<double>> total_g(pw.npw);
    std::vector<ModuleBase::Vector3<double>> gradient(pw.nrxx);
    std::vector<double> sigma(pw.nrxx, 0.0);
    std::vector<double> dedn(pw.nrxx, 0.0);
    std::vector<double> dedsigma(pw.nrxx, 0.0);
    std::vector<ModuleBase::Vector3<double>> field(pw.nrxx);
    std::vector<double> divergence(pw.nrxx, 0.0);
    pw.real2recip(total.data(), total_g.data());
    XC_Functional::grad_rho(total_g.data(), gradient.data(), &pw, pw.tpiba);
    for (int r = 0; r < pw.nrxx; ++r)
    {
        sigma[r] = gradient[r].x * gradient[r].x + gradient[r].y * gradient[r].y + gradient[r].z * gradient[r].z;
    }
    std::vector<double> direct_density = total;
    const double direct_energy = vdwxc_calculate(direct, direct_density.data(), sigma.data(), dedn.data(), dedsigma.data());
    for (int r = 0; r < pw.nrxx; ++r)
    {
        const double coefficient = 2.0 * dedsigma[r];
        field[r].x = coefficient * gradient[r].x;
        field[r].y = coefficient * gradient[r].y;
        field[r].z = coefficient * gradient[r].z;
    }
    XC_Functional::grad_dot(field.data(), divergence.data(), &pw, pw.tpiba);

    const double dv = pw.omega / pw.nxyz;
    double expected_vtxc = 0.0;
    EXPECT_NEAR(adapted.energy, 2.0 * direct_energy, 2.0e-12);
    for (int r = 0; r < pw.nrxx; ++r)
    {
        const double expected_potential = 2.0 * (dedn[r] - divergence[r]);
        expected_vtxc += dv * valence[r] * expected_potential;
        EXPECT_NEAR(adapted.potential[r], expected_potential, 2.0e-12);
    }
    EXPECT_NEAR(adapted.vtxc, expected_vtxc, 2.0e-12);
    vdwxc_finalize(&direct);
}
#endif

} // namespace
