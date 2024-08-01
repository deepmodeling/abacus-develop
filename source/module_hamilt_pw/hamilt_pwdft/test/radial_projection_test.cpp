#include "module_hamilt_pw/hamilt_pwdft/radial_projection.h"
#include <gtest/gtest.h>
#include <algorithm>
#include <numeric>
#include <fftw3.h>

#define DOUBLETHRESHOLD 1e-15
TEST(RadialProjectionTest, MaskfunctionGenerationTest)
{
    std::vector<double> mask;
    RadialProjection::_mask_func(mask);
    EXPECT_EQ(mask.size(), 201);
    EXPECT_EQ(mask[0], 1.0); // the rescaled value of the mask function, at 0, is 1
    EXPECT_NEAR(mask[200], 0.98138215E-05, 1e-10); // real space cut, at rc, is 0
}

TEST(RadialProjectionTest, ExampleTest)
{
    RadialProjection::RadialProjector rp;

    // use mask function as the example
    std::vector<double> mask;
    RadialProjection::_mask_func(mask);
    // suppose the r from 0 to 2.0 (inclusive) with 0.01 step
    std::vector<double> r(mask.size());
    std::iota(r.begin(), r.end(), 0);
    std::for_each(r.begin(), r.end(), [](double& x){x = x * 0.01;});

    // fill into the container
    std::vector<std::vector<double>> radials(1, mask);

    // suppose the angular momentum of function to transform is 0
    // but this is true because the mask function is really spherically
    // symmetric.
    std::vector<int> l(1, 0);

    // build the interpolation table
    rp._build_sbt_tab(r, radials, l, 201, 0.01); // build an interpolation table
    
    // only one q point: (0, 0, 0), is Gamma
    std::vector<ModuleBase::Vector3<double>> q(1);
    q[0] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);

    // so the output should have the same size as q
    std::vector<std::complex<double>> out(q.size());

    // transform, the last argument is 'r' for <G+k|p>.
    // for q = 0, the integration is actually performed on the function
    // itself, so it is actually the integration of the function, on the
    // grid from 0 to rc (2.0) with stepsize 0.01.
    rp.sbtft(q, out, 'r', 1.0, 1.0);
    // print, each 5 numbers in a row
    
    /**
     * The following Python code is used to generate the reference data
     * from scipy.integrate import simps
     * import numpy as np
     * 
     * r = np.arange(0, 2.01, 0.01)
     * mask = np.array(mask)
     * omega = 1
     * # Y00 is 1/sqrt(4*pi)
     * val = simps(mask*r**2, x=r) * np.sqrt(4*np.pi)
     * print(val/np.sqrt(omega))
     */
    std::complex<double> ref(0.6143605159327622, 0.0);
    EXPECT_NEAR(out[0].real(), ref.real(), DOUBLETHRESHOLD);
}

int main()
{
    testing::InitGoogleTest();
    return RUN_ALL_TESTS();
}