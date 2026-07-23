#include "gtest/gtest.h"
#include "xctest.h"
#include "../xc_functional.h"
#include "../exx_info.h"
#include "xc3_mock.h"
#include "source_base/matrix.h"

class XCTest_LaplacianAnalytical : public XCTest
{
  protected:
    ModulePW::PW_Basis rhopw;
    const int npw = 5;
    const int nrxx = 5;
    const int nmaxgr = 5;
    const double tpiba = 1.0;

    void SetUp() override
    {
        rhopw.nrxx = nrxx;
        rhopw.npw = npw;
        rhopw.nmaxgr = nmaxgr;
        rhopw.gcar = new ModuleBase::Vector3<double>[npw];
        for (int ig = 0; ig < npw; ig++)
            rhopw.gcar[ig] = ModuleBase::Vector3<double>(1.0, 1.0, 1.0);
    }

    void TearDown() override
    {
        delete[] rhopw.gcar;
    }
};

TEST_F(XCTest_LaplacianAnalytical, zero_input)
{
    std::vector<std::complex<double>> rhog(npw, 0.0);
    std::vector<double> lapl(nrxx, -1.0);
    XC_Functional::laplacian_rho(rhog.data(), lapl.data(), &rhopw, tpiba);
    for (int ir = 0; ir < nrxx; ir++)
        EXPECT_DOUBLE_EQ(lapl[ir], 0.0);
}

TEST_F(XCTest_LaplacianAnalytical, imaginary_rhog_constant)
{
    std::vector<std::complex<double>> rhog(npw);
    for (int ig = 0; ig < npw; ig++)
        rhog[ig] = std::complex<double>(0.0, 1.0);
    std::vector<double> lapl(nrxx);
    XC_Functional::laplacian_rho(rhog.data(), lapl.data(), &rhopw, tpiba);
    double g2 = 3.0;
    double expected = -g2 * tpiba * tpiba;
    for (int ir = 0; ir < nrxx; ir++)
        EXPECT_NEAR(lapl[ir], expected, 1e-14);
}

TEST_F(XCTest_LaplacianAnalytical, linearity)
{
    std::vector<std::complex<double>> rhog_a(npw);
    std::vector<std::complex<double>> rhog_b(npw);
    std::vector<std::complex<double>> rhog_sum(npw);
    for (int ig = 0; ig < npw; ig++)
    {
        rhog_a[ig] = std::complex<double>(ig, 2.0 * ig);
        rhog_b[ig] = std::complex<double>(3.0 * ig, ig);
        rhog_sum[ig] = rhog_a[ig] + rhog_b[ig];
    }

    std::vector<double> lapl_a(nrxx), lapl_b(nrxx), lapl_sum(nrxx);
    XC_Functional::laplacian_rho(rhog_a.data(), lapl_a.data(), &rhopw, tpiba);
    XC_Functional::laplacian_rho(rhog_b.data(), lapl_b.data(), &rhopw, tpiba);
    XC_Functional::laplacian_rho(rhog_sum.data(), lapl_sum.data(), &rhopw, tpiba);

    for (int ir = 0; ir < nrxx; ir++)
        EXPECT_NEAR(lapl_sum[ir], lapl_a[ir] + lapl_b[ir], 1e-14);
}

TEST_F(XCTest_LaplacianAnalytical, single_plane_wave)
{
    int n = 5;
    rhopw.nrxx = n;
    rhopw.npw = n;
    rhopw.nmaxgr = n;
    delete[] rhopw.gcar;
    rhopw.gcar = new ModuleBase::Vector3<double>[n];
    for (int ig = 0; ig < n; ig++)
        rhopw.gcar[ig] = ModuleBase::Vector3<double>(static_cast<double>(ig), 0.0, 0.0);

    std::vector<std::complex<double>> rhog(n, 0.0);
    rhog[0] = std::complex<double>(0.0, 1.0);
    std::vector<double> lapl(n);
    XC_Functional::laplacian_rho(rhog.data(), lapl.data(), &rhopw, tpiba);

    double g2 = 0.0;
    for (int i = 0; i < 3; i++)
        g2 += rhopw.gcar[0][i] * rhopw.gcar[0][i];
    double expected = -g2 * tpiba * tpiba;
    for (int ir = 0; ir < n; ir++)
        EXPECT_NEAR(lapl[ir], expected, 1e-14);
}
