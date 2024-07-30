#include "module_hamilt_pw/hamilt_pwdft/radial_projection.h"
#include <gtest/gtest.h>
#include <algorithm>
#include <numeric>
TEST(RadialProjectionTest, SBTFFTVector)
{
    // for a constant function, only the first element of the output array should be non-zero
    const int nr = 100;
    const std::vector<double> in(nr, 1.0); // constant function

    std::vector<double> r(nr); // r = 0, 0.01, 0.02, ..., 0.99
    std::iota(r.begin(), r.end(), 0);
    std::for_each(r.begin(), r.end(), [](double& x){x *= 0.01;});

    const int npw = 1;
    std::vector<ModuleBase::Vector3<double>> qs(npw, ModuleBase::Vector3<double>(1, 0, 0));
    RadialProjection::RadialProjector rp(3, qs, 1.0, 1.0);
    std::vector<std::complex<double>> out;
    rp.sbtfft(r, in, 2, out);
}

int main()
{
    testing::InitGoogleTest();
    return RUN_ALL_TESTS();
}