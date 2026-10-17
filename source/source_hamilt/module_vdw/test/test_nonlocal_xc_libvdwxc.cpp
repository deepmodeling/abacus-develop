#include "../nonlocal_xc_libvdwxc.h"

#include "source_basis/module_pw/pw_basis.h"
#include "source_hamilt/module_xc/xc_functional.h"

#include <cmath>
#include <complex>
#include <gtest/gtest.h>
#include <limits>
#include <type_traits>
#include <vector>

#ifdef __LIBVDWXC
// Direct library calls provide the API-parity reference. Include the C header
// after the C++ headers and drop the I macro that complex.h defines.
#include <complex.h>
#include <vdwxc.h>
#ifdef I
#undef I
#endif
#endif

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

    std::vector<double> perturbation() const
    {
        std::vector<double> values(pw.nrxx, 0.0);
        for (int x = 0; x < pw.nx; ++x)
            for (int y = 0; y < pw.ny; ++y)
                for (int z = 0; z < pw.nz; ++z)
                {
                    const int index = (x * pw.ny + y) * pw.nz + z;
                    const double ax = 2.0 * pi * x / pw.nx;
                    const double ay = 2.0 * pi * y / pw.ny;
                    values[index] = std::sin(ax + ay) + 0.3 * std::cos(2.0 * ay);
                }
        return values;
    }
};

#ifdef __LIBVDWXC
TEST(LibvdwxcOwnershipTest, BackendHandleCannotBeCopied)
{
    EXPECT_FALSE(std::is_copy_constructible<vdw::LibvdwxcNonlocalXC>::value);
    EXPECT_FALSE(std::is_copy_assignable<vdw::LibvdwxcNonlocalXC>::value);
}

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
                        // libvdwxc stores each C00,C10,C20 group as one
                        // lattice vector; ABACUS stores lattice vectors in rows.
                        a.e11 * pw.lat0, a.e12 * pw.lat0, a.e13 * pw.lat0,
                        a.e21 * pw.lat0, a.e22 * pw.lat0, a.e23 * pw.lat0,
                        a.e31 * pw.lat0, a.e32 * pw.lat0, a.e33 * pw.lat0);
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

TEST_F(LibvdwxcAdapterTest, RejectsLeftHandedCell)
{
    ModulePW::PW_Basis reflected("cpu", "double");
    const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -2.0);
    reflected.initgrids(10.0, lattice, 6, 8, 12);
    reflected.initparameters(false, 4.0, 1, true);
    reflected.setuptransform();
    reflected.collect_local_pw();
    EXPECT_EXIT(vdw::LibvdwxcNonlocalXC(reflected, 6.3, 0.0093), testing::ExitedWithCode(1), "");
}

TEST_F(LibvdwxcAdapterTest, RejectsInvalidParametersWithStandardDiagnostic)
{
    EXPECT_EXIT(vdw::LibvdwxcNonlocalXC(pw, 0.0, 0.0093), testing::ExitedWithCode(1), "");
}

TEST_F(LibvdwxcAdapterTest, RejectsInvalidDensityWithStandardDiagnostic)
{
    vdw::LibvdwxcNonlocalXC adapter(pw, 6.3, 0.0093);
    const std::vector<double> total = density(0.004, 0.006);
    std::vector<double> invalid = total;
    invalid.pop_back();
    EXPECT_EXIT(adapter.evaluate(pw, invalid, total), testing::ExitedWithCode(1), "");
    invalid = total;
    invalid[0] = std::numeric_limits<double>::quiet_NaN();
    EXPECT_EXIT(adapter.evaluate(pw, invalid, total), testing::ExitedWithCode(1), "");
}

TEST_F(LibvdwxcAdapterTest, CellMappingPreservesIndependentReciprocalGeometry)
{
    // Both cells have volume sqrt(3)*lat0^3 and |G100|=4*pi/(sqrt(3)*lat0).
    // A density depending only on fractional x must therefore give identical
    // isotropic nonlocal energy; transposing the oblique cell breaks this equality.
    ModulePW::PW_Basis orthogonal("cpu", "double");
    const ModuleBase::Matrix3 lattice(std::sqrt(0.75), 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 2.0);
    orthogonal.initgrids(10.0, lattice, pw.nx, pw.ny, pw.nz);
    orthogonal.initparameters(false, 4.0, 1, true);
    orthogonal.setuptransform();
    orthogonal.collect_local_pw();

    const double amplitude = 0.003;
    std::vector<double> total(pw.nrxx);
    for (int x = 0; x < pw.nx; ++x)
        for (int y = 0; y < pw.ny; ++y)
            for (int z = 0; z < pw.nz; ++z)
            {
                const int index = (x * pw.ny + y) * pw.nz + z;
                const double phase = 2.0 * pi * x / pw.nx;
                total[index] = 0.01 + amplitude * std::cos(phase);
            }

    // The hand-derived gradient of cos(2*pi*s1) uses b1=(1,1/sqrt(3),0).
    std::vector<std::complex<double>> density_g(pw.npw);
    std::vector<ModuleBase::Vector3<double>> gradient(pw.nrxx);
    pw.real2recip(total.data(), density_g.data());
    XC_Functional::grad_rho(density_g.data(), gradient.data(), &pw, pw.tpiba);
    for (int x = 0; x < pw.nx; ++x)
        for (int y = 0; y < pw.ny; ++y)
            for (int z = 0; z < pw.nz; ++z)
            {
                const int index = (x * pw.ny + y) * pw.nz + z;
                const double phase = 2.0 * pi * x / pw.nx;
                const double derivative = -amplitude * (2.0 * pi / pw.lat0) * std::sin(phase);
                EXPECT_NEAR(gradient[index].x, derivative, 2.0e-12);
                EXPECT_NEAR(gradient[index].y, derivative / std::sqrt(3.0), 2.0e-12);
                EXPECT_NEAR(gradient[index].z, 0.0, 2.0e-12);
            }

    vdw::LibvdwxcNonlocalXC oblique_adapter(pw, 6.3, 0.0093);
    vdw::LibvdwxcNonlocalXC orthogonal_adapter(orthogonal, 6.3, 0.0093);
    const vdw::NonlocalXCResult oblique_result = oblique_adapter.evaluate(pw, total, total);
    const vdw::NonlocalXCResult orthogonal_result = orthogonal_adapter.evaluate(orthogonal, total, total);
    EXPECT_NEAR(pw.omega, orthogonal.omega, 2.0e-12);
    ASSERT_NEAR(oblique_result.energy, orthogonal_result.energy, 2.0e-12);
    EXPECT_NEAR(oblique_result.vtxc, orthogonal_result.vtxc, 2.0e-12);
    for (int r = 0; r < pw.nrxx; ++r)
    {
        EXPECT_NEAR(oblique_result.potential[r], orthogonal_result.potential[r], 2.0e-12);
    }
}

TEST_F(LibvdwxcAdapterTest, PotentialMatchesFiniteDifferenceOfEnergy)
{
    const double b = 6.3;
    const double c = 0.0093;
    const double step = 1.0e-5;
    const std::vector<double> total = density(0.004, 0.006);
    const std::vector<double> valence = total;
    const std::vector<double> direction = perturbation();
    vdw::LibvdwxcNonlocalXC adapter(pw, b, c);
    const vdw::NonlocalXCResult reference = adapter.evaluate(pw, total, valence);
    std::vector<double> plus = total;
    std::vector<double> minus = total;
    for (int r = 0; r < pw.nrxx; ++r)
    {
        plus[r] += step * direction[r];
        minus[r] -= step * direction[r];
    }
    const vdw::NonlocalXCResult forward = adapter.evaluate(pw, plus, valence);
    const vdw::NonlocalXCResult backward = adapter.evaluate(pw, minus, valence);
    const double finite_difference = (forward.energy - backward.energy) / (2.0 * step);
    const double dv = pw.omega / pw.nxyz;
    double potential_contraction = 0.0;
    for (int r = 0; r < pw.nrxx; ++r)
    {
        potential_contraction += dv * reference.potential[r] * direction[r];
    }
    EXPECT_NEAR(finite_difference, potential_contraction, 2.0e-6);
}
#endif

} // namespace
