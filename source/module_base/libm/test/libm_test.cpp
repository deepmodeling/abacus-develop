#include <vector>
#include <complex>
#include <cmath>
#include <random>
#include <float.h>

#include "gtest/gtest.h"
#include "../libm.h"

TEST(base_libm, sincos)
{
    int len = 50000;
    std::vector<double> da(len);
    std::vector<double> ds1(len*2);
    std::vector<double> ds2(len*2);

    std::uniform_real_distribution<double> rnd(-50000, 50000);
    std::default_random_engine eng;

    for (int i = 0; i < len; ++i) {
        da[i] = rnd(eng);
    }

    for (int i = 0; i < len; ++i) {
        sincos(da[i], &ds1[i * 2 + 0], &ds1[i * 2 + 1]);
        ModuleBase::libm::sincos(da[i], &ds2[i * 2 + 0], &ds2[i * 2 + 1]);
    }

    for (int i = 0; i < len * 2; i++) {
        EXPECT_DOUBLE_EQ(ds1[i], ds2[i]);
    }
}

TEST(base_libm, exp)
{
    int len = 50000;
    std::vector<double> da(len);
    std::vector<double> ds1(len);
    std::vector<double> ds2(len);

    std::uniform_real_distribution<double> rnd(-1000, 1000);
    std::default_random_engine eng;

    for (int i = 0; i < len; ++i) {
        da[i] = rnd(eng);
    }

    for (int i = 0; i < len; ++i) {
        ds1[i] = std::exp(da[i]);
        ds2[i] = ModuleBase::libm::exp(da[i]);
    }

    for (int i = 0; i < len; i++) {
        EXPECT_DOUBLE_EQ(ds1[i], ds2[i]);
    }
}

TEST(base_libm, cexp)
{
    int len = 50000;
    std::vector<std::complex<double>> da(len);
    std::vector<std::complex<double>> ds1(len);
    std::vector<std::complex<double>> ds2(len);

    std::uniform_real_distribution<double> rnd(-1000, 1000);
    std::default_random_engine eng;

    for (int i = 0; i < len; ++i) {
        da[i] = std::complex<double>(rnd(eng), rnd(eng));
    }

    for (int i = 0; i < len; ++i) {
        ds1[i] = std::exp(da[i]);
        ds2[i] = ModuleBase::libm::exp(da[i]);
    }

    for (int i = 0; i < len; i++) {
        EXPECT_DOUBLE_EQ(ds1[i].real(), ds2[i].real());
        EXPECT_DOUBLE_EQ(ds1[i].imag(), ds2[i].imag());
    }
}
