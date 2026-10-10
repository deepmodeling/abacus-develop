#include "../xc_functional.h"
#include "../libxc_abacus.h"
#include "gtest/gtest.h"
#include "xctest.h"
#include <array>
#include <cmath>
#include <iostream>
#include <iomanip>

namespace ModuleBase
{
    void WARNING_QUIT(const std::string &file,const std::string &description) {exit(1);}
    void TITLE(const std::string &class_function_name,bool disable){};
    void TITLE(const std::string &class_name,const std::string &function_name,bool disable){};
}

class XCTest_SCANL_Laplacian : public XCTest
{
    protected:
        double e_base, v1_base, v2_base, v3_base, vlapl_base;
        double e_modified, v1_modified, v2_modified, v3_modified, vlapl_modified;
        double e_scaled, v1_scaled, v2_scaled, v3_scaled, vlapl_scaled;

        void SetUp()
        {
            XC_Functional::set_xc_type("MGGA_X_SCANL+MGGA_C_SCANL");

            const double rho = 0.17E+01;
            const double grho = 0.81E-11;
            const double tau = 0.02403590412;
            const double lapl_base = 0.15E+01;
            double hybrid_alpha = 0.0;
            double hse_omega = 0.0;

            XC_Functional_Libxc::tau_xc(
                XC_Functional::get_func_id(),
                rho, grho, lapl_base, tau,
                e_base, v1_base, v2_base, v3_base, vlapl_base, hybrid_alpha, hse_omega
            );

            XC_Functional_Libxc::tau_xc(
                XC_Functional::get_func_id(),
                rho, grho, lapl_base + 1.0, tau,
                e_modified, v1_modified, v2_modified, v3_modified, vlapl_modified, hybrid_alpha, hse_omega
            );

            XC_Functional_Libxc::tau_xc(
                XC_Functional::get_func_id(),
                rho, grho, 2.0 * lapl_base, tau,
                e_scaled, v1_scaled, v2_scaled, v3_scaled, vlapl_scaled, hybrid_alpha, hse_omega
            );
        }
};

TEST_F(XCTest_SCANL_Laplacian, laplacian_affects_energy)
{
    EXPECT_NE(e_base, e_modified);
    EXPECT_NE(e_base, e_scaled);

    std::cout << std::scientific << std::setprecision(15);
    std::cout << "\n=== SCAN Laplacian Sensitivity Test ===" << std::endl;
    std::cout << "Base Laplacian:    " << 0.15E+01 << std::endl;
    std::cout << "  E_xc =           " << e_base << std::endl;
    std::cout << "  dE/dlapl ≈       " << (e_modified - e_base) / 1.0 << std::endl;
    std::cout << "Modified Laplacian:" << 0.15E+01 + 1.0 << std::endl;
    std::cout << "  E_xc =           " << e_modified << std::endl;
    std::cout << "  Delta E =        " << e_modified - e_base << std::endl;
    std::cout << "Scaled Laplacian:  " << 2.0 * 0.15E+01 << std::endl;
    std::cout << "  E_xc =           " << e_scaled << std::endl;
    std::cout << "  Delta E =        " << e_scaled - e_base << std::endl;
    std::cout << "========================================" << std::endl;
}

// Verify direct libxc call returns finite physical values.
// Reference values are libxc-version-dependent; asserting against hardcoded
// numbers would break on other libxc versions (e.g. 7.1.2). Instead we check
// that the values are finite and the wrapper matches the direct libxc call.
TEST(XC_ScanL_Reference, MGGA_X_SCANL_direct_libxc)
{
    const double rho   = 35.536521214608185;
    const double sigma = 1.149382202334535e+05;
    const double lapl  = 8.411855859277239e+02;
    const double tau   = 1.389887953757970e+02;

    std::vector<int> func_id = {XC_MGGA_X_SCANL};
    std::vector<xc_func_type> funcs = XC_Functional_Libxc::init_func(func_id, XC_UNPOLARIZED, 0.0, 0.0);

    double exc = 0.0, vrho = 0.0, vsigma = 0.0, vlapl = 0.0, vtau = 0.0;
    xc_mgga_exc_vxc(&funcs[0], 1, &rho, &sigma, &lapl, &tau,
                    &exc, &vrho, &vsigma, &vlapl, &vtau);

    // SCAN-L exchange: negative energy, non-vanishing vlapl, zero vtau
    EXPECT_LT(exc, 0.0);
    EXPECT_TRUE(std::isfinite(exc));
    EXPECT_TRUE(std::isfinite(vrho));
    EXPECT_TRUE(std::isfinite(vsigma));
    EXPECT_NE(vlapl, 0.0);
    EXPECT_DOUBLE_EQ(vtau, 0.0);

    XC_Functional_Libxc::finish_func(funcs);
}

TEST(XC_ScanL_Reference, MGGA_C_SCANL_direct_libxc)
{
    const double rho   = 35.536521214608185;
    const double sigma = 1.149382202334535e+05;
    const double lapl  = 8.411855859277239e+02;
    const double tau   = 1.389887953757970e+02;

    std::vector<int> func_id = {XC_MGGA_C_SCANL};
    std::vector<xc_func_type> funcs = XC_Functional_Libxc::init_func(func_id, XC_UNPOLARIZED, 0.0, 0.0);

    double exc = 0.0, vrho = 0.0, vsigma = 0.0, vlapl = 0.0, vtau = 0.0;
    xc_mgga_exc_vxc(&funcs[0], 1, &rho, &sigma, &lapl, &tau,
                    &exc, &vrho, &vsigma, &vlapl, &vtau);

    // SCAN-L correlation: negative energy, non-vanishing vlapl, zero vtau
    EXPECT_LT(exc, 0.0);
    EXPECT_TRUE(std::isfinite(exc));
    EXPECT_TRUE(std::isfinite(vrho));
    EXPECT_TRUE(std::isfinite(vsigma));
    EXPECT_NE(vlapl, 0.0);
    EXPECT_DOUBLE_EQ(vtau, 0.0);

    XC_Functional_Libxc::finish_func(funcs);
}

// Verify tau_xc wrapper produces the same results as the direct libxc call.
// This validates the ABACUS wrapper without depending on a specific libxc
// version's absolute reference values.
TEST(XC_ScanL_Reference, tau_xc_wrapper_libxc)
{
    XC_Functional::set_xc_type("MGGA_X_SCANL+MGGA_C_SCANL");

    const double rho   = 35.536521214608185;
    const double grho  = 1.149382202334535e+05;
    const double lapl  = 8.411855859277239e+02;
    const double tau   = 1.389887953757970e+02;

    // Direct libxc calls for exchange and correlation
    std::vector<int> func_id = XC_Functional::get_func_id();
    std::vector<xc_func_type> funcs = XC_Functional_Libxc::init_func(func_id, XC_UNPOLARIZED, 0.0, 0.0);
    double exc_x = 0.0, vrho_x = 0.0, vsigma_x = 0.0, vlapl_x = 0.0, vtau_x = 0.0;
    xc_mgga_exc_vxc(&funcs[0], 1, &rho, &grho, &lapl, &tau,
                    &exc_x, &vrho_x, &vsigma_x, &vlapl_x, &vtau_x);
    double exc_c = 0.0, vrho_c = 0.0, vsigma_c = 0.0, vlapl_c = 0.0, vtau_c = 0.0;
    xc_mgga_exc_vxc(&funcs[1], 1, &rho, &grho, &lapl, &tau,
                    &exc_c, &vrho_c, &vsigma_c, &vlapl_c, &vtau_c);
    XC_Functional_Libxc::finish_func(funcs);

    double sxc_ref = (exc_x + exc_c) * rho;

    double sxc, v1xc, v2xc, v3xc, vlaplxc;
    double hybrid_alpha = 0.0;
    double hse_omega = 0.0;
    XC_Functional_Libxc::tau_xc(XC_Functional::get_func_id(), rho, grho, lapl, tau,
                                sxc, v1xc, v2xc, v3xc, vlaplxc, hybrid_alpha, hse_omega);

    EXPECT_NEAR(sxc, sxc_ref, std::abs(sxc_ref) * 1.0e-6);
    EXPECT_NE(v1xc, 0.0);
    EXPECT_NE(vlaplxc, 0.0);
}

#ifdef __EXX
// SCAN0 scales the semilocal SCAN exchange by (1 - hybrid_alpha). The
// wrappers must use the hybrid_alpha they are given, not the global
// XC_Functional value, so the two are deliberately set to different values.
class XCTest_SCAN0_Alpha : public testing::Test
{
  protected:
    const std::vector<int> x_only = {XC_MGGA_X_SCAN};
    const double global_alpha = 0.25;
    const double hybrid_alpha = 0.4;
    const double hse_omega = 0.0;

    void SetUp() override
    {
        XC_Functional::set_xc_type("SCAN0");
        XC_Functional::set_hybrid_alpha(global_alpha);
    }
};

TEST_F(XCTest_SCAN0_Alpha, tau_xc)
{
    const double no_alpha = 0.0;
    const double scale = 1.0 - hybrid_alpha;
    const double rho = 0.5;
    const double grho = 0.01;
    const double lapl = 0.0;
    const double tau = 0.1;

    std::array<double, 5> full = {0.0, 0.0, 0.0, 0.0, 0.0};
    XC_Functional_Libxc::tau_xc(x_only, rho, grho, lapl, tau,
                                full[0], full[1], full[2], full[3], full[4], no_alpha, hse_omega);

    std::array<double, 5> scaled = {0.0, 0.0, 0.0, 0.0, 0.0};
    XC_Functional_Libxc::tau_xc(x_only, rho, grho, lapl, tau,
                                scaled[0], scaled[1], scaled[2], scaled[3], scaled[4], hybrid_alpha, hse_omega);

    std::vector<xc_func_type> funcs = XC_Functional_Libxc::init_func(x_only, XC_UNPOLARIZED, hybrid_alpha, hse_omega);
    std::array<double, 5> scaled_funcs = {0.0, 0.0, 0.0, 0.0, 0.0};
    XC_Functional_Libxc::tau_xc(funcs, rho, grho, lapl, tau,
                                scaled_funcs[0], scaled_funcs[1], scaled_funcs[2], scaled_funcs[3], scaled_funcs[4],
                                hybrid_alpha);
    XC_Functional_Libxc::finish_func(funcs);

    EXPECT_NE(full[0], 0.0);
    for (int i = 0; i < 5; i++)
    {
        EXPECT_NEAR(scaled[i], scale * full[i], 1.0e-12);
        EXPECT_NEAR(scaled_funcs[i], scale * full[i], 1.0e-12);
    }
}

TEST_F(XCTest_SCAN0_Alpha, tau_xc_spin)
{
    const double no_alpha = 0.0;
    const double scale = 1.0 - hybrid_alpha;
    const double rhoup = 0.3;
    const double rhodw = 0.2;
    const ModuleBase::Vector3<double> gdr1(0.05, 0.02, 0.01);
    const ModuleBase::Vector3<double> gdr2(0.03, 0.01, 0.02);
    const double laplup = 0.0;
    const double lapldw = 0.0;
    const double tauup = 0.06;
    const double taudw = 0.04;

    std::array<double, 10> full = {};
    XC_Functional_Libxc::tau_xc_spin(x_only, rhoup, rhodw, gdr1, gdr2, laplup, lapldw, tauup, taudw,
                                     full[0], full[1], full[2], full[3], full[4], full[5],
                                     full[6], full[7], full[8], full[9], no_alpha, hse_omega);

    std::array<double, 10> scaled = {};
    XC_Functional_Libxc::tau_xc_spin(x_only, rhoup, rhodw, gdr1, gdr2, laplup, lapldw, tauup, taudw,
                                     scaled[0], scaled[1], scaled[2], scaled[3], scaled[4], scaled[5],
                                     scaled[6], scaled[7], scaled[8], scaled[9], hybrid_alpha, hse_omega);

    EXPECT_NE(full[0], 0.0);
    for (int i = 0; i < 10; i++)
    {
        EXPECT_NEAR(scaled[i], scale * full[i], 1.0e-12);
    }
}
#endif
