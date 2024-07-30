#include "module_hamilt_pw/hamilt_pwdft/radial_projection.h"
#include <gtest/gtest.h>
#include <algorithm>
#include <numeric>

TEST(RadialProjectionTest, MaskfunctionGenerationTest)
{
    std::vector<double> mask;
    RadialProjection::_mask_func(mask);
    EXPECT_EQ(mask.size(), 201);
}

TEST(RadialProjectionTest, ExampleTest)
{
    RadialProjection::RadialProjector rp;
    // use mask function as the example
    std::vector<double> mask;
    RadialProjection::_mask_func(mask);
    std::vector<double> r(mask.size());
    std::iota(r.begin(), r.end(), 0);
    std::for_each(r.begin(), r.end(), [](double& x){x = x * 0.01;});

    std::vector<std::vector<double>> radials(1, mask);
    std::vector<int> l(1, 0);

    rp._build_sbt_tab(r, radials, l, 200, 0.005); // build an interpolation table
    std::vector<double> result;

    std::vector<ModuleBase::Vector3<double>> q(1);
    q[0] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
    std::vector<std::complex<double>> out;
    rp.sbtft(q, out, 1.0, 1.0, 'r');

}

int main()
{
    testing::InitGoogleTest();
    return RUN_ALL_TESTS();
}