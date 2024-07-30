#include "module_hamilt_pw/hamilt_pwdft/radial_projection.h"
#include <gtest/gtest.h>
#include <algorithm>
#include <numeric>
TEST(RadialProjectionTest, SBTFFTVector)
{
    // for a constant function, only the first element of the output array should be non-zero
    const int l = 2;
    const int m = 0;
    const int nr = 100;
    const std::vector<double> in(nr, 1.0); // constant function
    const double omega = 1;
    const double tpiba = 1;

    std::vector<double> r(nr); // r = 0.01, 0.02, 0.03, ...
    std::iota(r.begin(), r.end(), 1);
    std::transform(r.begin(), r.end(), r.begin(), [](double x) { return x * 0.01; });
    std::complex<double> out;
    const ModuleBase::Vector3<double> q{0.0, 0.0, 0.0};

    RadialProjection::_sbtfft(in, r, l, m, q, omega, tpiba, out);
    EXPECT_GT(std::abs(out), 1e-10);
}

int main()
{
    testing::InitGoogleTest();
    return RUN_ALL_TESTS();
}