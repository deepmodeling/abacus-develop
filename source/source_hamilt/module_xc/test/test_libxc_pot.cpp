#include "../libxc_abacus.h"
#include "../xc_functional.h"
#include "source_estate/module_charge/charge.h"

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <tuple>
#include <vector>

namespace
{
// This target uses the one-grid-point FFT mocks defined in test_xc5.cpp.
// All storage here is owned by the fixture, including the arrays borrowed by
// the mocked Charge and PW_Basis classes (whose destructors are empty).
struct MetaPoint
{
    ModulePW::PW_Basis basis;
    Charge charge;
    std::array<ModuleBase::Vector3<double>, 1> reciprocal;
    std::array<double, 2> density;
    std::array<double, 2> kinetic;
    std::array<double*, 2> density_rows;
    std::array<double*, 2> kinetic_rows;
    std::array<double, 1> core;

    MetaPoint()
    {
        reciprocal[0] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
        density = {{0.2, 0.1}};
        kinetic = {{0.6, 0.4}};
        core = {{0.0}};
        density_rows = {{&density[0], &density[1]}};
        kinetic_rows = {{&kinetic[0], &kinetic[1]}};
        basis.nrxx = 1;
        basis.npw = 1;
        basis.nxyz = 1;
        basis.nmaxgr = 1;
        basis.gcar = reciprocal.data();
        charge.rhopw = &basis;
        charge.rho = density_rows.data();
        charge.kin_r = kinetic_rows.data();
        charge.rho_core = core.data();
    }
};

typedef std::tuple<double, double, ModuleBase::matrix, ModuleBase::matrix> MetaResult;

MetaResult evaluate_meta(MetaPoint& point, const int nspin)
{
    const std::vector<int> functionals = {XC_MGGA_X_SCAN, XC_MGGA_C_SCAN};
    const double omega = 1.0;
    const double tpiba = 1.0;
    const double hybrid_alpha = 0.0;
    const double hse_omega = 0.0;
    return XC_Functional_Libxc::v_xc_meta(functionals, point.basis.nrxx, omega, tpiba,
                                         &point.charge, nspin, hybrid_alpha, hse_omega);
}

void check_uniform_limit(const int nspin)
{
    MetaPoint point;
    const MetaResult result = evaluate_meta(point, nspin);
    const int polarization = nspin == 1 ? XC_UNPOLARIZED : XC_POLARIZED;
    const std::array<double, 3> sigma = {{0.0, 0.0, 0.0}};
    const std::array<double, 2> laplacian = {{0.0, 0.0}};
    const std::array<double, 2> tau = {{0.3, 0.2}}; // ABACUS stores 2 * tau.
    double expected_energy = 0.0;
    std::array<double, 2> expected_potential = {{0.0, 0.0}};
    std::array<double, 2> expected_tau_potential = {{0.0, 0.0}};
    const std::array<int, 2> ids = {{XC_MGGA_X_SCAN, XC_MGGA_C_SCAN}};
    for (const int id : ids)
    {
        xc_func_type functional;
        ASSERT_EQ(xc_func_init(&functional, id, polarization), 0);
        double energy = 0.0;
        std::array<double, 2> vrho = {{0.0, 0.0}};
        std::array<double, 3> vsigma = {{0.0, 0.0, 0.0}};
        std::array<double, 2> vlapl = {{0.0, 0.0}};
        std::array<double, 2> vtau = {{0.0, 0.0}};
        xc_mgga_exc_vxc(&functional, 1, point.density.data(), sigma.data(), laplacian.data(),
                        tau.data(), &energy, vrho.data(), vsigma.data(), vlapl.data(), vtau.data());
        xc_func_end(&functional);
        for (int spin = 0; spin < nspin; ++spin)
        {
            expected_energy += 2.0 * energy * point.density[spin];
            expected_potential[spin] += 2.0 * vrho[spin];
            expected_tau_potential[spin] += vtau[spin];
        }
    }
    EXPECT_LT(expected_energy, -1e-3);
    EXPECT_NEAR(std::get<0>(result), expected_energy, 1e-12);
    double expected_vtxc = 0.0;
    for (int spin = 0; spin < nspin; ++spin)
    {
        EXPECT_NEAR(std::get<2>(result)(spin, 0), expected_potential[spin], 1e-12);
        EXPECT_NEAR(std::get<3>(result)(spin, 0), expected_tau_potential[spin], 1e-12);
        expected_vtxc += expected_potential[spin] * point.density[spin];
    }
    EXPECT_NEAR(std::get<1>(result), expected_vtxc, 1e-12);
}

class LibxcMetaPotentialTest : public testing::Test
{
  protected:
    void SetUp() override
    {
        XC_Functional::set_xc_type("SCAN");
    }
};
} // namespace

TEST_F(LibxcMetaPotentialTest, UniformUnpolarizedMatchesLibxc)
{
    check_uniform_limit(1);
}

TEST_F(LibxcMetaPotentialTest, UniformPolarizedMatchesLibxc)
{
    check_uniform_limit(2);
}

TEST_F(LibxcMetaPotentialTest, LocalTermsRemainContinuousAcrossGradientCutoff)
{
    for (const int nspin : {1, 2})
    {
        MetaPoint point;
        const MetaResult zero = evaluate_meta(point, nspin);
        point.reciprocal[0].x = 1e-13;
        const MetaResult below = evaluate_meta(point, nspin);
        point.reciprocal[0].x = 1e-10;
        const MetaResult above = evaluate_meta(point, nspin);
        EXPECT_NEAR(std::get<0>(zero), std::get<0>(above), 1e-12);
        EXPECT_NEAR(std::get<0>(below), std::get<0>(above), 1e-12);
        for (int spin = 0; spin < nspin; ++spin)
        {
            EXPECT_NEAR(std::get<2>(zero)(spin, 0), std::get<2>(above)(spin, 0), 1e-12);
            EXPECT_NEAR(std::get<3>(zero)(spin, 0), std::get<3>(above)(spin, 0), 1e-12);
        }
    }
}

TEST_F(LibxcMetaPotentialTest, ExistingDensityAndKineticCutoffsArePreserved)
{
    for (const int nspin : {1, 2})
    {
        MetaPoint low_density;
        low_density.density = {{1e-12, 1e-12}};
        const MetaResult density_result = evaluate_meta(low_density, nspin);
        EXPECT_DOUBLE_EQ(std::get<0>(density_result), 0.0);
        MetaPoint low_kinetic;
        low_kinetic.kinetic = {{1e-12, 1e-12}};
        const MetaResult kinetic_result = evaluate_meta(low_kinetic, nspin);
        EXPECT_DOUBLE_EQ(std::get<0>(kinetic_result), 0.0);
    }
}
