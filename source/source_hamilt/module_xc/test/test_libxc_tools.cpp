#include "../xc_functional.h"
#include "../libxc_abacus.h"
#include "gtest/gtest.h"
#include "xctest.h"
#include "../exx_info.h"
#include "xc3_mock.h"
#include "source_base/matrix.h"
#include "source_cell/cal_ux.h"
#include "../../../source_base/parallel_reduce.h"

/************************************************
 *  unit test of XC_Functional_Libxc::cal_sgn_vxc
 *  (libxc_tools.cpp)
 *
 *  cal_sgn_vxc returns two threshold masks following the
 *  convention of Quantum ESPRESSO's libxc interface
 *  (XClib/xc_wrapper_gga.f90):
 *  - the first mask applies to exc and vrho, zeroed only
 *    where the density falls below rho_threshold_vrho;
 *  - the second mask applies only to the vsigma (gradient)
 *    term, zeroed where the density falls below
 *    rho_threshold_vsigma or sqrt(|sigma|) falls below
 *    grho_threshold_vsigma.
 *  For non-GGA families the vsigma mask must stay 1.
 ***********************************************/

namespace
{
// same values as the call site in v_xc_libxc (libxc_pot.cpp)
constexpr double rho_threshold_vrho = 1E-10;
constexpr double rho_threshold_vsigma = 1E-6;
constexpr double grho_threshold_vsigma = 1E-10;
} // namespace

class LibxcToolsSgnVxcTest : public XCTest
{
  protected:
    xc_func_type gga_func;
    xc_func_type lda_func;

    void SetUp() override
    {
        xc_func_init(&gga_func, XC_GGA_X_PBE, XC_UNPOLARIZED);
        xc_func_init(&lda_func, XC_LDA_X, XC_UNPOLARIZED);
    }

    void TearDown() override
    {
        xc_func_end(&gga_func);
        xc_func_end(&lda_func);
    }
};

// nspin = 1, GGA functional: the vrho mask survives in the shell
// rho_threshold_vrho < rho <= rho_threshold_vsigma, while the vsigma
// mask is already zeroed there
TEST_F(LibxcToolsSgnVxcTest, Nspin1GgaDensityTiers)
{
    const std::size_t nrxx = 3;
    const std::vector<double> rho = {1.0, 1E-8, 1E-12};
    const std::vector<double> sigma = {0.1, 0.1, 0.1};

    const std::pair<std::vector<double>, std::vector<double>> sgn = XC_Functional_Libxc::cal_sgn_vxc(
        rho_threshold_vrho, rho_threshold_vsigma, grho_threshold_vsigma, gga_func, 1, nrxx, rho, sigma);

    EXPECT_DOUBLE_EQ(sgn.first[0], 1.0);
    EXPECT_DOUBLE_EQ(sgn.second[0], 1.0);
    // inside the shell: vrho kept, vsigma suppressed
    EXPECT_DOUBLE_EQ(sgn.first[1], 1.0);
    EXPECT_DOUBLE_EQ(sgn.second[1], 0.0);
    // below rho_threshold_vrho: both suppressed
    EXPECT_DOUBLE_EQ(sgn.first[2], 0.0);
    EXPECT_DOUBLE_EQ(sgn.second[2], 0.0);
}

// nspin = 1, GGA functional: a small gradient alone suppresses
// vsigma even where the density is large
TEST_F(LibxcToolsSgnVxcTest, Nspin1GgaGradientTrigger)
{
    const std::size_t nrxx = 2;
    const std::vector<double> rho = {1.0, 1.0};
    const std::vector<double> sigma = {0.1, 1E-22}; // sqrt(1E-22) = 1E-11 < grho_threshold_vsigma

    const std::pair<std::vector<double>, std::vector<double>> sgn = XC_Functional_Libxc::cal_sgn_vxc(
        rho_threshold_vrho, rho_threshold_vsigma, grho_threshold_vsigma, gga_func, 1, nrxx, rho, sigma);

    EXPECT_DOUBLE_EQ(sgn.first[0], 1.0);
    EXPECT_DOUBLE_EQ(sgn.second[0], 1.0);
    EXPECT_DOUBLE_EQ(sgn.first[1], 1.0);
    EXPECT_DOUBLE_EQ(sgn.second[1], 0.0);
}

// nspin = 1, LDA functional: no vsigma term exists, so the vsigma
// mask must stay 1 regardless of density and sigma; the vrho mask
// still follows rho_threshold_vrho
TEST_F(LibxcToolsSgnVxcTest, Nspin1LdaKeepsVsigmaMask)
{
    const std::size_t nrxx = 2;
    const std::vector<double> rho = {1E-8, 1E-12};
    const std::vector<double> sigma = {0.0, 0.0};

    const std::pair<std::vector<double>, std::vector<double>> sgn = XC_Functional_Libxc::cal_sgn_vxc(
        rho_threshold_vrho, rho_threshold_vsigma, grho_threshold_vsigma, lda_func, 1, nrxx, rho, sigma);

    EXPECT_DOUBLE_EQ(sgn.first[0], 1.0);
    EXPECT_DOUBLE_EQ(sgn.second[0], 1.0);
    EXPECT_DOUBLE_EQ(sgn.first[1], 0.0);
    EXPECT_DOUBLE_EQ(sgn.second[1], 1.0);
}

// nspin = 2, GGA functional: both spin channels are masked jointly
// when either spin density fails a threshold
TEST_F(LibxcToolsSgnVxcTest, Nspin2JointMasking)
{
    const std::size_t nrxx = 3;
    // interleaved spin densities: {up0, dw0, up1, dw1, up2, dw2}
    const std::vector<double> rho = {1E-8, 0.5, 1E-12, 0.5, 0.5, 0.5};
    // interleaved sigma: {uu0, ud0, dd0, uu1, ud1, dd1, ...}
    const std::vector<double> sigma = {0.1, 0.0, 0.1, 0.1, 0.0, 0.1, 0.1, 0.0, 1E-22};

    const std::pair<std::vector<double>, std::vector<double>> sgn = XC_Functional_Libxc::cal_sgn_vxc(
        rho_threshold_vrho, rho_threshold_vsigma, grho_threshold_vsigma, gga_func, 2, nrxx, rho, sigma);

    // ir0: rho_up inside the shell -> vsigma suppressed in both channels,
    // vrho kept in both channels
    EXPECT_DOUBLE_EQ(sgn.first[0], 1.0);
    EXPECT_DOUBLE_EQ(sgn.first[1], 1.0);
    EXPECT_DOUBLE_EQ(sgn.second[0], 0.0);
    EXPECT_DOUBLE_EQ(sgn.second[1], 0.0);
    // ir1: rho_up below rho_threshold_vrho -> both masks suppressed
    // in both channels
    EXPECT_DOUBLE_EQ(sgn.first[2], 0.0);
    EXPECT_DOUBLE_EQ(sgn.first[3], 0.0);
    EXPECT_DOUBLE_EQ(sgn.second[2], 0.0);
    EXPECT_DOUBLE_EQ(sgn.second[3], 0.0);
    // ir2: densities fine but the down-down gradient fails ->
    // only vsigma suppressed, in both channels
    EXPECT_DOUBLE_EQ(sgn.first[4], 1.0);
    EXPECT_DOUBLE_EQ(sgn.first[5], 1.0);
    EXPECT_DOUBLE_EQ(sgn.second[4], 0.0);
    EXPECT_DOUBLE_EQ(sgn.second[5], 0.0);
}

// nspin = 2, GGA functional: the up-down cross component of sigma
// does not trigger the vsigma mask by itself
TEST_F(LibxcToolsSgnVxcTest, Nspin2CrossSigmaIgnored)
{
    const std::size_t nrxx = 1;
    const std::vector<double> rho = {0.5, 0.5};
    const std::vector<double> sigma = {0.1, 1E-30, 0.1};

    const std::pair<std::vector<double>, std::vector<double>> sgn = XC_Functional_Libxc::cal_sgn_vxc(
        rho_threshold_vrho, rho_threshold_vsigma, grho_threshold_vsigma, gga_func, 2, nrxx, rho, sigma);

    EXPECT_DOUBLE_EQ(sgn.first[0], 1.0);
    EXPECT_DOUBLE_EQ(sgn.first[1], 1.0);
    EXPECT_DOUBLE_EQ(sgn.second[0], 1.0);
    EXPECT_DOUBLE_EQ(sgn.second[1], 1.0);
}
