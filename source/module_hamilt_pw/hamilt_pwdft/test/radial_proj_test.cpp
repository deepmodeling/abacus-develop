#include "module_hamilt_pw/hamilt_pwdft/radial_proj.h"
#include <gtest/gtest.h>
#include <algorithm>
#include <numeric>
#include <fftw3.h>
#include <random>

#define DOUBLETHRESHOLD 1e-15

TEST(RadialProjectionTest, BuildBackwardMapTest)
{
    const std::vector<std::vector<int>> it2iproj = {
        {0, 1},
        {2, 3, 4, 5},
        {6, 7, 8}
    };
    const std::vector<int> iproj2l = {0, 1, 0, 1, 2, 3, 0, 1, 2};
    /**
     * a backward map will index irow to (it, iproj, m). For the above example,
     * it is actually the first atom type has two projectors, the first projector
     * has l = 0, the second has l = 1. The former will have only one channel, 
     * the second will have 3. The second atom type has 4 projectors...
     * 
     * Therefore, it is expected there are 1+3+1+3+5+7+1+3+5 = 29 rows in total.
     */
    std::vector<int> itref(4, 0); // four 0, sixteen 1, nine 2
    itref.insert(itref.end(), 16, 1);
    itref.insert(itref.end(), 9, 2);
    std::vector<int> iprojref(1, 0);
    iprojref.insert(iprojref.end(), 3, 1);
    iprojref.insert(iprojref.end(), 1, 0);
    iprojref.insert(iprojref.end(), 3, 1);
    iprojref.insert(iprojref.end(), 5, 2);
    iprojref.insert(iprojref.end(), 7, 3);
    iprojref.insert(iprojref.end(), 1, 0);
    iprojref.insert(iprojref.end(), 3, 1);
    iprojref.insert(iprojref.end(), 5, 2);
    std::vector<int> mref(1, 0);
    mref.insert(mref.end(), {-1, 0, 1});
    mref.insert(mref.end(), {0});
    mref.insert(mref.end(), {-1, 0, 1});
    mref.insert(mref.end(), {-2, -1, 0, 1, 2});
    mref.insert(mref.end(), {-3, -2, -1, 0, 1, 2, 3});
    mref.insert(mref.end(), {0});
    mref.insert(mref.end(), {-1, 0, 1});
    mref.insert(mref.end(), {-2, -1, 0, 1, 2});

    std::vector<int> irow2it, irow2iproj, irow2m;
    RadialProjection::RadialProjector::_build_backward_map(it2iproj, iproj2l, irow2it, irow2iproj, irow2m);
    for(int i = 0; i < 29; ++i)
    {
        EXPECT_EQ(irow2it[i], itref[i]);
        EXPECT_EQ(irow2iproj[i], iprojref[i]);
        EXPECT_EQ(irow2m[i], mref[i]);
    }
}

TEST(RadialProjectionTest, BuildForwardMapTest)
{
    // forward means from (it, ia, iproj, m) to irow
    const std::vector<std::vector<int>> it2ia = {
        {0, 1},
        {0, 1, 3},
        {12}
    };
    const std::vector<std::vector<int>> it2iproj = {
        {0, 1},
        {2, 3, 4, 5},
        {6, 7, 8}
    };
    const std::vector<int> iproj2l = {0, 1, 0, 1, 2, 3, 0, 1, 2};
    std::map<std::tuple<int, int, int, int>, int> forward_map;
    RadialProjection::RadialProjector::_build_forward_map(it2ia, it2iproj, iproj2l, forward_map);
    // this is a terribly large map, for type 0, there are 2 atoms, each atom has two projectors
    // with angular momentum as 0 and 1, respectively. For type 1, there are 3 atoms, each has
    // 4 projectors with angular momentum from 0 to 3...
    const int nrow = 2*4 + 3*16 + 1*9;
    EXPECT_EQ(forward_map.size(), nrow);
    //                                  it  ia iproj m
    // it = 0
    // ia = 0
    // iproj = 0: l = 0
    EXPECT_EQ(forward_map[std::make_tuple(0, 0, 0, 0)], 0);
    // iproj = 1: l = 1: m = -1, 0, 1
    EXPECT_EQ(forward_map[std::make_tuple(0, 0, 1, -1)], 1);
    EXPECT_EQ(forward_map[std::make_tuple(0, 0, 1, 0)], 2);
    EXPECT_EQ(forward_map[std::make_tuple(0, 0, 1, 1)], 3);
    // ia = 1
    // iproj = 0: l = 0
    EXPECT_EQ(forward_map[std::make_tuple(0, 1, 0, 0)], 4);
    // iproj = 1: l = 1: m = -1, 0, 1
    EXPECT_EQ(forward_map[std::make_tuple(0, 1, 1, -1)], 5);
    EXPECT_EQ(forward_map[std::make_tuple(0, 1, 1, 0)], 6);
    EXPECT_EQ(forward_map[std::make_tuple(0, 1, 1, 1)], 7);
    // it = 1
    // ia = 0
    // iproj = 0: l = 0
    EXPECT_EQ(forward_map[std::make_tuple(1, 0, 0, 0)], 8);
    // iproj = 1: l = 1: m = -1, 0, 1
    EXPECT_EQ(forward_map[std::make_tuple(1, 0, 1, -1)], 9);
    EXPECT_EQ(forward_map[std::make_tuple(1, 0, 1, 0)], 10);
    EXPECT_EQ(forward_map[std::make_tuple(1, 0, 1, 1)], 11);
    // iproj = 2: l = 2: m = -2, -1, 0, 1, 2
    EXPECT_EQ(forward_map[std::make_tuple(1, 0, 2, -2)], 12);
    EXPECT_EQ(forward_map[std::make_tuple(1, 0, 2, -1)], 13);
    EXPECT_EQ(forward_map[std::make_tuple(1, 0, 2, 0)], 14);
    EXPECT_EQ(forward_map[std::make_tuple(1, 0, 2, 1)], 15);
    EXPECT_EQ(forward_map[std::make_tuple(1, 0, 2, 2)], 16);
    // iproj = 3: l = 3: m = -3, -2, -1, 0, 1, 2, 3
    EXPECT_EQ(forward_map[std::make_tuple(1, 0, 3, -3)], 17);
    EXPECT_EQ(forward_map[std::make_tuple(1, 0, 3, -2)], 18);
    EXPECT_EQ(forward_map[std::make_tuple(1, 0, 3, -1)], 19);
    EXPECT_EQ(forward_map[std::make_tuple(1, 0, 3, 0)], 20);
    EXPECT_EQ(forward_map[std::make_tuple(1, 0, 3, 1)], 21);
    EXPECT_EQ(forward_map[std::make_tuple(1, 0, 3, 2)], 22);
    EXPECT_EQ(forward_map[std::make_tuple(1, 0, 3, 3)], 23);
    // ia = 1...
}

TEST(RadialProjectionTest, MaskfunctionGenerationTest)
{
    std::vector<double> mask;
    RadialProjection::_mask_func(mask);
    EXPECT_EQ(mask.size(), 201);
    EXPECT_EQ(mask[0], 1.0); // the rescaled value of the mask function, at 0, is 1
    EXPECT_NEAR(mask[200], 0.98138215E-05, 1e-10); // real space cut, at rc, is 0
}

TEST(RadialProjectionTest, BuildSbtTabCorrectnessTest)
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

TEST(RadialProjectionTest, BuildSbtTabStabilityTest)
{
    // still use mask function but scale with different Gaussian function
    std::vector<double> mask;
    RadialProjection::_mask_func(mask);
    // suppose the r from 0 to 2.0 (inclusive) with 0.01 step
    std::vector<double> r(mask.size());
    std::iota(r.begin(), r.end(), 0);
    std::for_each(r.begin(), r.end(), [](double& x){x = x * 0.01;});

    std::vector<double> aug1(mask.size());
    std::vector<double> aug2(mask.size());
    std::vector<double> aug3(mask.size());
    const double sigma1 = 0.1;
    const double sigma2 = 1.0;
    const double sigma3 = 10.0;
    std::transform(r.begin(), r.end(), aug1.begin(), [sigma1](double x){return std::exp(-x*x/(2*sigma1*sigma1));});
    std::transform(r.begin(), r.end(), aug2.begin(), [sigma2](double x){return std::exp(-x*x/(2*sigma2*sigma2));});
    std::transform(r.begin(), r.end(), aug3.begin(), [sigma3](double x){return std::exp(-x*x/(2*sigma3*sigma3));});

    std::vector<double> mask1(mask.size());
    std::vector<double> mask2(mask.size());
    std::vector<double> mask3(mask.size());
    std::transform(mask.begin(), mask.end(), aug1.begin(), mask1.begin(), std::multiplies<double>());
    std::transform(mask.begin(), mask.end(), aug2.begin(), mask2.begin(), std::multiplies<double>());
    std::transform(mask.begin(), mask.end(), aug3.begin(), mask3.begin(), std::multiplies<double>());

    // then prepare the q vectors, in total 50, each with 3 components among which each is random
    // number ranging from 0 to 2, use the random seed as 1
    const int npw = 50;
    std::vector<ModuleBase::Vector3<double>> q(npw);
    std::mt19937 gen(1);
    std::uniform_real_distribution<double> dis(0.0, 2.0);
    std::for_each(q.begin(), q.end(), [&gen, &dis](ModuleBase::Vector3<double>& x){
        x[0] = dis(gen);
        x[1] = dis(gen);
        x[2] = dis(gen);
    });

    // fill into the container
    std::vector<std::vector<double>> radials(4); // do not remember mask0 = mask
    radials[0] = mask;
    radials[1] = mask1;
    radials[2] = mask2;
    radials[3] = mask3;

    // angular momentum set as 0, 1, 2, 3
    std::vector<int> l(4);
    std::iota(l.begin(), l.end(), 0);

    RadialProjection::RadialProjector rp;
    rp._build_sbt_tab(r, radials, l, 401, 0.01); // build an interpolation table

    // then perform the transform
    std::vector<std::complex<double>> out;
    rp.sbtft(q, out, 'r', 1.0, 1.0);

    // there are 50 q points, and radials have l as 0, 1, 2, 3, which means the
    // number of lm components in total is 1 + 3 + 5 + 7 = 16
    // therefore there should be 50 * 16 = 800 complex numbers in the output
    EXPECT_EQ(out.size(), npw*16);
    // check if they are listed in expected sequence, say the first should be
    // Fourier transform of the first function, mask, l = 0, m = 0
    RadialProjection::RadialProjector rp1;
    rp1._build_sbt_tab(r, 
                       std::vector<std::vector<double>>(1, mask), 
                       std::vector<int>(1, 0), 401, 0.01);
    std::vector<std::complex<double>> out1;
    rp1.sbtft(q, out1, 'r', 1.0, 1.0);
    for(int iq = 0; iq < 50; ++iq) // 50 q-points
    {
        EXPECT_NEAR(out[iq].real(), out1[iq].real(), DOUBLETHRESHOLD);
        EXPECT_NEAR(out[iq].imag(), out1[iq].imag(), DOUBLETHRESHOLD);
    }
    RadialProjection::RadialProjector rp2;
    rp2._build_sbt_tab(r, 
                       std::vector<std::vector<double>>(1, mask1), 
                       std::vector<int>(1, 1), 401, 0.01);
    std::vector<std::complex<double>> out2;
    rp2.sbtft(q, out2, 'r', 1.0, 1.0);
    for(int iq = 0; iq < 50*(2*1+1); ++iq) // 50 q-points
    {
        EXPECT_NEAR(out[iq+50].real(), out2[iq].real(), DOUBLETHRESHOLD);
        EXPECT_NEAR(out[iq+50].imag(), out2[iq].imag(), DOUBLETHRESHOLD);
    }
    RadialProjection::RadialProjector rp3;
    rp3._build_sbt_tab(r, 
                       std::vector<std::vector<double>>(1, mask2), 
                       std::vector<int>(1, 2), 401, 0.01);
    std::vector<std::complex<double>> out3;
    rp3.sbtft(q, out3, 'r', 1.0, 1.0);
    for(int iq = 0; iq < 50*(2*2+1); ++iq) // 50 q-points
    {
        EXPECT_NEAR(out[iq+50+50*(2*1+1)].real(), out3[iq].real(), DOUBLETHRESHOLD);
        EXPECT_NEAR(out[iq+50+50*(2*1+1)].imag(), out3[iq].imag(), DOUBLETHRESHOLD);
    }
    RadialProjection::RadialProjector rp4;
    rp4._build_sbt_tab(r, 
                       std::vector<std::vector<double>>(1, mask3), 
                       std::vector<int>(1, 3), 401, 0.01);
    std::vector<std::complex<double>> out4;
    rp4.sbtft(q, out4, 'r', 1.0, 1.0);
    for(int iq = 0; iq < 50*(2*3+1); ++iq) // 50 q-points
    {
        EXPECT_NEAR(out[iq+50+50*(2*1+1)+50*(2*2+1)].real(), out4[iq].real(), DOUBLETHRESHOLD);
        EXPECT_NEAR(out[iq+50+50*(2*1+1)+50*(2*2+1)].imag(), out4[iq].imag(), DOUBLETHRESHOLD);
    }
}

int main()
{
    testing::InitGoogleTest();
    return RUN_ALL_TESTS();
}