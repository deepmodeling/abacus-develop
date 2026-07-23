#include "../xc_functional.h"
#include "../libxc_abacus.h"
#include "gtest/gtest.h"
#include "xctest.h"
#include "../exx_info.h"
#include <cmath>
#include <iostream>
#include <iomanip>

namespace ModuleBase
{
    void WARNING_QUIT(const std::string &file,const std::string &description) {exit(1);}
    void TITLE(const std::string &class_function_name,bool disable){};
    void TITLE(const std::string &class_name,const std::string &function_name,bool disable){};
}

namespace GlobalV
{
    std::string BASIS_TYPE = "";
    bool CAL_STRESS = false;
    int CAL_FORCE = 0;
    int NSPIN = 1;
}

namespace GlobalC
{
	Exx_Info exx_info;
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
    std::cout << "\n=== SCAN-L Laplacian Sensitivity Test ===" << std::endl;
    std::cout << "Base Laplacian:    " << 0.15E+01 << std::endl;
    std::cout << "  E_xc =           " << e_base << std::endl;
    std::cout << "  dE/dlapl ~=      " << (e_modified - e_base) / 1.0 << std::endl;
    std::cout << "Modified Laplacian:" << 0.15E+01 + 1.0 << std::endl;
    std::cout << "  E_xc =           " << e_modified << std::endl;
    std::cout << "  Delta E =        " << e_modified - e_base << std::endl;
    std::cout << "Scaled Laplacian:  " << 2.0 * 0.15E+01 << std::endl;
    std::cout << "  E_xc =           " << e_scaled << std::endl;
    std::cout << "  Delta E =        " << e_scaled - e_base << std::endl;
    std::cout << "=========================================" << std::endl;
}

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

    EXPECT_NEAR(exc,   -2.569101943337681e+00, 5.0e-04);
    EXPECT_NEAR(vrho,  -2.912514021641506e+00, 1.0e-03);
    EXPECT_NEAR(vsigma, -8.289754258579972e-05, 1.0e-12);
    EXPECT_NEAR(vlapl,  3.971745124876924e-03, 1.0e-12);
    EXPECT_EQ(vtau, 0.0);

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

    EXPECT_NEAR(exc,   -5.250261832501450e-02, 1.0e-06);
    EXPECT_NEAR(vrho,  -1.323740339176898e-01, 1.0e-05);
    EXPECT_NEAR(vsigma, 1.138021340988220e-05, 1.0e-10);
    EXPECT_NEAR(vlapl, -3.867564402989276e-04, 1.0e-10);
    EXPECT_EQ(vtau, 0.0);

    XC_Functional_Libxc::finish_func(funcs);
}

TEST(XC_ScanL_Reference, tau_xc_wrapper_libxc)
{
    XC_Functional::set_xc_type("MGGA_X_SCANL+MGGA_C_SCANL");

    const double rho   = 35.536521214608185;
    const double grho  = 1.149382202334535e+05;
    const double lapl  = 8.411855859277239e+02;
    const double tau   = 1.389887953757970e+02;

    double sxc, v1xc, v2xc, v3xc, vlaplxc;
    double hybrid_alpha = 0.0;
    double hse_omega = 0.0;
    XC_Functional_Libxc::tau_xc(XC_Functional::get_func_id(), rho, grho, lapl, tau,
                                sxc, v1xc, v2xc, v3xc, vlaplxc, hybrid_alpha, hse_omega);

    double exc_x_ref = -2.569101943337681e+00;
    double exc_c_ref = -5.250261832501450e-02;
    double sxc_ref   = (exc_x_ref + exc_c_ref) * rho;

    EXPECT_NEAR(sxc, sxc_ref, 0.05);
    EXPECT_NE(v1xc, 0.0);
    EXPECT_NE(vlaplxc, 0.0);
}
