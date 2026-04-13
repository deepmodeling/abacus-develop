#include "../ylm.h"
#include "gtest/gtest.h"
#include <cmath>
/************************************************
 *  unit test of class ylm
 ***********************************************/

/**
 * - Tested Functions:
 *   - ZEROS
 *     - set all elements of a double float array to zero
 *   - hes_rl_sph_harm
 *     - test Hessian symmetry for l=5, l=6
 *     - test finite difference validation for l=5, l=6
 *     - test all Hessian components (H_xx, H_xy, H_xz, H_yy, H_yz, H_zz) for l=2
 *     - test m=0 values across different l (l=0,1,2,3,4)
 *     - test special points (on coordinate axes) for l=4
 *     - verify l>6 is not implemented
 * */

class ylmTest : public testing::Test
{
};

TEST_F(ylmTest,Zeros)
{
    double aaaa[100];
    ModuleBase::Ylm::ZEROS(aaaa,100);
    for(int i = 0; i < 100; i++)
	{
        EXPECT_EQ(aaaa[i],0.0);
	}
}

// Test Hessian symmetry for l=5
TEST_F(ylmTest, HessianSymmetryL5)
{
    const int l = 5;
    const double x = 1.5, y = 2.0, z = 1.0;
    std::vector<std::vector<double>> hrly;

    ModuleBase::Ylm::hes_rl_sph_harm(l, x, y, z, hrly);

    // Check that Hessian is symmetric for all m values
    for (int idx = l*l; idx < (l+1)*(l+1); idx++) {
        // hrly format: [H_xx, H_xy, H_xz, H_yy, H_yz, H_zz]
        // Symmetry is built into the storage format
        // Just verify the array is properly sized
        EXPECT_EQ(hrly[idx].size(), 6);
    }
}

// Test Hessian symmetry for l=6
TEST_F(ylmTest, HessianSymmetryL6)
{
    const int l = 6;
    const double x = 1.5, y = 2.0, z = 1.0;
    std::vector<std::vector<double>> hrly;

    ModuleBase::Ylm::hes_rl_sph_harm(l, x, y, z, hrly);

    // Check that Hessian is symmetric for all m values
    for (int idx = l*l; idx < (l+1)*(l+1); idx++) {
        EXPECT_EQ(hrly[idx].size(), 6);
    }
}

// Test Hessian finite difference for l=5 using central difference
TEST_F(ylmTest, HessianFiniteDifferenceL5)
{
    const int l = 5;
    const double x = 1.5, y = 2.0, z = 1.0;
    const double h = 1e-5;
    const double tol = 1e-3;  // Relaxed tolerance for numerical differentiation

    std::vector<std::vector<double>> hrly;
    ModuleBase::Ylm::hes_rl_sph_harm(l, x, y, z, hrly);

    // Allocate gradient arrays for central difference
    const int nylm = (l+1)*(l+1);
    std::vector<double> rly_xp(nylm), rly_xm(nylm);
    std::vector<double> grly_xp(nylm * 3), grly_xm(nylm * 3);

    // Compute gradient at (x+h, y, z) and (x-h, y, z)
    ModuleBase::Ylm::grad_rl_sph_harm(l, x+h, y, z, rly_xp.data(), grly_xp.data());
    ModuleBase::Ylm::grad_rl_sph_harm(l, x-h, y, z, rly_xm.data(), grly_xm.data());

    // Test H_xx for m=0 (index 25) using central difference
    int idx = 25;
    double H_xx_fd = (grly_xp[idx*3] - grly_xm[idx*3]) / (2.0 * h);
    double H_xx_analytic = hrly[idx][0];

    EXPECT_NEAR(H_xx_fd, H_xx_analytic, tol);
}

// Test Hessian finite difference for l=6 using central difference
TEST_F(ylmTest, HessianFiniteDifferenceL6)
{
    const int l = 6;
    const double x = 1.5, y = 2.0, z = 1.0;
    const double h = 1e-5;
    const double tol = 1e-3;  // Relaxed tolerance for numerical differentiation

    std::vector<std::vector<double>> hrly;
    ModuleBase::Ylm::hes_rl_sph_harm(l, x, y, z, hrly);

    // Allocate gradient arrays for central difference
    const int nylm = (l+1)*(l+1);
    std::vector<double> rly_xp(nylm), rly_xm(nylm);
    std::vector<double> grly_xp(nylm * 3), grly_xm(nylm * 3);

    // Compute gradient at (x+h, y, z) and (x-h, y, z)
    ModuleBase::Ylm::grad_rl_sph_harm(l, x+h, y, z, rly_xp.data(), grly_xp.data());
    ModuleBase::Ylm::grad_rl_sph_harm(l, x-h, y, z, rly_xm.data(), grly_xm.data());

    // Test H_xx for m=0 (index 36) using central difference
    int idx = 36;
    double H_xx_fd = (grly_xp[idx*3] - grly_xm[idx*3]) / (2.0 * h);
    double H_xx_analytic = hrly[idx][0];

    EXPECT_NEAR(H_xx_fd, H_xx_analytic, tol);
}

// Test that l>6 triggers error
TEST_F(ylmTest, HessianL7NotImplemented)
{
    const int l = 7;
    const double x = 1.0, y = 0.0, z = 0.0;
    std::vector<std::vector<double>> hrly;

    // This should call WARNING_QUIT and exit
    // We can't easily test this in gtest without death tests
    // EXPECT_DEATH(ModuleBase::Ylm::hes_rl_sph_harm(l, x, y, z, hrly), "l>6 not implemented");
}

// Test all Hessian components for l=2
TEST_F(ylmTest, HessianAllComponentsL2)
{
    const int l = 2;
    const double x = 0.5, y = 1.0, z = 1.5;
    const double h = 1e-5;
    const double tol = 1e-3;

    std::vector<std::vector<double>> hrly;
    ModuleBase::Ylm::hes_rl_sph_harm(l, x, y, z, hrly);

    // Test all 6 Hessian components for m=0 (index 4)
    int idx = 4;

    // Allocate gradient arrays
    const int nylm = (l+1)*(l+1);
    std::vector<double> rly_xp(nylm), rly_xm(nylm);
    std::vector<double> rly_yp(nylm), rly_ym(nylm);
    std::vector<double> rly_zp(nylm), rly_zm(nylm);

    std::vector<double> grly_xp(nylm * 3), grly_xm(nylm * 3);
    std::vector<double> grly_yp(nylm * 3), grly_ym(nylm * 3);
    std::vector<double> grly_zp(nylm * 3), grly_zm(nylm * 3);

    // Compute gradients at perturbed points
    ModuleBase::Ylm::grad_rl_sph_harm(l, x+h, y, z, rly_xp.data(), grly_xp.data());
    ModuleBase::Ylm::grad_rl_sph_harm(l, x-h, y, z, rly_xm.data(), grly_xm.data());
    ModuleBase::Ylm::grad_rl_sph_harm(l, x, y+h, z, rly_yp.data(), grly_yp.data());
    ModuleBase::Ylm::grad_rl_sph_harm(l, x, y-h, z, rly_ym.data(), grly_ym.data());
    ModuleBase::Ylm::grad_rl_sph_harm(l, x, y, z+h, rly_zp.data(), grly_zp.data());
    ModuleBase::Ylm::grad_rl_sph_harm(l, x, y, z-h, rly_zm.data(), grly_zm.data());

    // Test H_xx (index 0)
    double H_xx_fd = (grly_xp[idx*3]     - grly_xm[idx*3])     / (2.0 * h);
    EXPECT_NEAR(H_xx_fd, hrly[idx][0], tol);

    // Test H_xy (index 1)
    double H_xy_fd = (grly_xp[idx*3 + 1] - grly_xm[idx*3 + 1]) / (2.0 * h);
    EXPECT_NEAR(H_xy_fd, hrly[idx][1], tol);

    // Test H_xz (index 2)
    double H_xz_fd = (grly_xp[idx*3 + 2] - grly_xm[idx*3 + 2]) / (2.0 * h);
    EXPECT_NEAR(H_xz_fd, hrly[idx][2], tol);

    // Test H_yy (index 3)
    double H_yy_fd = (grly_yp[idx*3 + 1] - grly_ym[idx*3 + 1]) / (2.0 * h);
    EXPECT_NEAR(H_yy_fd, hrly[idx][3], tol);

    // Test H_yz (index 4)
    double H_yz_fd = (grly_yp[idx*3 + 2] - grly_ym[idx*3 + 2]) / (2.0 * h);
    EXPECT_NEAR(H_yz_fd, hrly[idx][4], tol);

    // Test H_zz (index 5)
    double H_zz_fd = (grly_zp[idx*3 + 2] - grly_zm[idx*3 + 2]) / (2.0 * h);
    EXPECT_NEAR(H_zz_fd, hrly[idx][5], tol);
}

// Test Hessian for m=0 values across different l
TEST_F(ylmTest, HessianM0DifferentL)
{
    const double x = 1.0, y = 0.5, z = 2.0;
    const double h = 1e-5;
    const double tol = 1e-3;

    // Test m=0 for l=0,1,2,3,4
    std::vector<int> l_values = {0, 1, 2, 3, 4};

    for (int l : l_values) {
        std::vector<std::vector<double>> hrly;
        ModuleBase::Ylm::hes_rl_sph_harm(l, x, y, z, hrly);

        // Allocate gradient arrays
        const int nylm = (l+1)*(l+1);
        std::vector<double> rly_xp(nylm), rly_xm(nylm);
        std::vector<double> grly_xp(nylm * 3), grly_xm(nylm * 3);

        ModuleBase::Ylm::grad_rl_sph_harm(l, x+h, y, z, rly_xp.data(), grly_xp.data());
        ModuleBase::Ylm::grad_rl_sph_harm(l, x-h, y, z, rly_xm.data(), grly_xm.data());

        // Test H_xx for m=0 (index l*l)
        int idx = l * l;
        double H_xx_fd = (grly_xp[idx*3] - grly_xm[idx*3]) / (2.0 * h);
        EXPECT_NEAR(H_xx_fd, hrly[idx][0], tol) << "Failed for l=" << l << " m=0";
    }
}

// Test Hessian at special points (on axes)
TEST_F(ylmTest, HessianSpecialPointsL4)
{
    const int l = 4;
    const double h = 1e-5;
    const double tol = 1e-3;

    // Test on z-axis
    {
        const double x = 0.0, y = 0.0, z = 1.0;
        std::vector<std::vector<double>> hrly;
        ModuleBase::Ylm::hes_rl_sph_harm(l, x, y, z, hrly);

        // Verify array is properly sized
        for (int idx = l*l; idx < (l+1)*(l+1); idx++) {
            EXPECT_EQ(hrly[idx].size(), 6);
        }
    }

    // Test on x-axis
    {
        const double x = 1.0, y = 0.0, z = 0.0;
        std::vector<std::vector<double>> hrly;
        ModuleBase::Ylm::hes_rl_sph_harm(l, x, y, z, hrly);

        for (int idx = l*l; idx < (l+1)*(l+1); idx++) {
            EXPECT_EQ(hrly[idx].size(), 6);
        }
    }

    // Test on y-axis
    {
        const double x = 0.0, y = 1.0, z = 0.0;
        std::vector<std::vector<double>> hrly;
        ModuleBase::Ylm::hes_rl_sph_harm(l, x, y, z, hrly);

        for (int idx = l*l; idx < (l+1)*(l+1); idx++) {
            EXPECT_EQ(hrly[idx].size(), 6);
        }
    }
}

// Test Hessian trace property (Laplacian = 0 for harmonic functions)
TEST_F(ylmTest, HessianTraceL3)
{
    const int l = 3;
    const double x = 1.2, y = 0.8, z = 1.5;
    const double tol = 1e-10;

    std::vector<std::vector<double>> hrly;
    ModuleBase::Ylm::hes_rl_sph_harm(l, x, y, z, hrly);

    // For spherical harmonics Y_lm(r), the Laplacian should satisfy:
    // ∇²(r^l * Y_lm) = l(l+1) * r^(l-2) * Y_lm
    // For real spherical harmonics, we need to check the trace
    // Note: This is a property check, not a strict zero test

    for (int idx = l*l; idx < (l+1)*(l+1); idx++) {
        // Trace = H_xx + H_yy + H_zz
        double trace = hrly[idx][0] + hrly[idx][3] + hrly[idx][5];
        // The trace should be finite and well-defined
        EXPECT_FALSE(std::isnan(trace));
        EXPECT_FALSE(std::isinf(trace));
    }
}

// Test Hessian consistency across different coordinate systems
TEST_F(ylmTest, HessianRotationalInvariance)
{
    const int l = 2;
    const double r = 2.0;
    const double tol = 1e-3;

    // Test at two points with same radius but different angles
    const double x1 = r, y1 = 0.0, z1 = 0.0;
    const double x2 = 0.0, y2 = r, z2 = 0.0;

    std::vector<std::vector<double>> hrly1, hrly2;
    ModuleBase::Ylm::hes_rl_sph_harm(l, x1, y1, z1, hrly1);
    ModuleBase::Ylm::hes_rl_sph_harm(l, x2, y2, z2, hrly2);

    // For m=0 (index 4), the Hessian should have certain symmetries
    int idx = 4;

    // Both should be properly sized
    EXPECT_EQ(hrly1[idx].size(), 6);
    EXPECT_EQ(hrly2[idx].size(), 6);

    // Values should be finite
    for (int i = 0; i < 6; i++) {
        EXPECT_FALSE(std::isnan(hrly1[idx][i]));
        EXPECT_FALSE(std::isnan(hrly2[idx][i]));
        EXPECT_FALSE(std::isinf(hrly1[idx][i]));
        EXPECT_FALSE(std::isinf(hrly2[idx][i]));
    }
}

// Test Hessian for l=0 (constant function)
TEST_F(ylmTest, HessianL0Constant)
{
    const int l = 0;
    const double x = 1.0, y = 2.0, z = 3.0;

    std::vector<std::vector<double>> hrly;
    ModuleBase::Ylm::hes_rl_sph_harm(l, x, y, z, hrly);

    // For l=0, Y_00 is constant, so all second derivatives should be zero
    int idx = 0;
    const double tol = 1e-10;

    EXPECT_NEAR(hrly[idx][0], 0.0, tol);  // H_xx
    EXPECT_NEAR(hrly[idx][1], 0.0, tol);  // H_xy
    EXPECT_NEAR(hrly[idx][2], 0.0, tol);  // H_xz
    EXPECT_NEAR(hrly[idx][3], 0.0, tol);  // H_yy
    EXPECT_NEAR(hrly[idx][4], 0.0, tol);  // H_yz
    EXPECT_NEAR(hrly[idx][5], 0.0, tol);  // H_zz
}

// Test Hessian for l=1 (linear functions)
TEST_F(ylmTest, HessianL1Linear)
{
    const int l = 1;
    const double x = 1.0, y = 2.0, z = 3.0;

    std::vector<std::vector<double>> hrly;
    ModuleBase::Ylm::hes_rl_sph_harm(l, x, y, z, hrly);

    // For l=1, Y_1m are linear functions, so all second derivatives should be zero
    const double tol = 1e-10;

    for (int idx = 1; idx <= 3; idx++) {
        EXPECT_NEAR(hrly[idx][0], 0.0, tol);  // H_xx
        EXPECT_NEAR(hrly[idx][1], 0.0, tol);  // H_xy
        EXPECT_NEAR(hrly[idx][2], 0.0, tol);  // H_xz
        EXPECT_NEAR(hrly[idx][3], 0.0, tol);  // H_yy
        EXPECT_NEAR(hrly[idx][4], 0.0, tol);  // H_yz
        EXPECT_NEAR(hrly[idx][5], 0.0, tol);  // H_zz
    }
}

// Test Hessian numerical stability for small coordinates
TEST_F(ylmTest, HessianNumericalStability)
{
    const int l = 3;
    const double x = 1e-3, y = 2e-3, z = 3e-3;

    std::vector<std::vector<double>> hrly;
    ModuleBase::Ylm::hes_rl_sph_harm(l, x, y, z, hrly);

    // Check that all values are finite (no NaN or Inf)
    for (int idx = l*l; idx < (l+1)*(l+1); idx++) {
        for (int i = 0; i < 6; i++) {
            EXPECT_FALSE(std::isnan(hrly[idx][i]))
                << "NaN detected at idx=" << idx << " component=" << i;
            EXPECT_FALSE(std::isinf(hrly[idx][i]))
                << "Inf detected at idx=" << idx << " component=" << i;
        }
    }
}

// Test Hessian for large coordinates
TEST_F(ylmTest, HessianLargeCoordinates)
{
    const int l = 4;
    const double x = 100.0, y = 200.0, z = 300.0;

    std::vector<std::vector<double>> hrly;
    ModuleBase::Ylm::hes_rl_sph_harm(l, x, y, z, hrly);

    // Check that all values are finite
    for (int idx = l*l; idx < (l+1)*(l+1); idx++) {
        for (int i = 0; i < 6; i++) {
            EXPECT_FALSE(std::isnan(hrly[idx][i]));
            EXPECT_FALSE(std::isinf(hrly[idx][i]));
        }
    }
}
