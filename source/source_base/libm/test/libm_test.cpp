#include <vector>
#include <complex>
#include <cmath>
#include <random>
#include <float.h>

#include "gtest/gtest.h"
#include "../libm.h"

#define MY_EXPECT_DOUBLE_EQ(ds1, ds2) \
do {                                                                            \
    if (std::isnan(ds1) && std::isnan(ds2) &&                                   \
        std::signbit(ds1) == std::signbit(ds2))                                 \
        EXPECT_EQ(1, 1);                                                        \
    else EXPECT_DOUBLE_EQ(ds1, ds2);                                            \
} while (0)

#define MY_EXPECT_COMPLEX_DOUBLE_EQ(ds1, ds2)                                   \
    MY_EXPECT_DOUBLE_EQ(ds1.real(), ds2.real());                                \
    MY_EXPECT_DOUBLE_EQ(ds1.imag(), ds2.imag());                                \

TEST(base_libm, sincos_random)
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
        MY_EXPECT_DOUBLE_EQ(ds1[i], ds2[i]);
    }
}

TEST(base_libm, sincos_spec)
{
    
    double da[] = {
        -0.0, +0.0,
        -INFINITY, +INFINITY,
        -NAN, +NAN,
        -DBL_MIN, +DBL_MIN,
        -DBL_MAX, +DBL_MAX,
        -DBL_MIN * DBL_MIN, +DBL_MIN * DBL_MIN,
        -DBL_MAX * DBL_MAX, +DBL_MAX * DBL_MAX,
        -DBL_MAX / M_PI};

    const int len = sizeof(da) / sizeof(double);
    std::vector<double> ds1(len*2);
    std::vector<double> ds2(len*2);

    for (int i = 0; i < len; ++i) {
        sincos(da[i], &ds1[i * 2 + 0], &ds1[i * 2 + 1]);
        ModuleBase::libm::sincos(da[i], &ds2[i * 2 + 0], &ds2[i * 2 + 1]);
    }

    for (int i = 0; i < len * 2; i++) {
        MY_EXPECT_DOUBLE_EQ(ds1[i], ds2[i]);
    }
}

#define SINGLE_PARAM_FLOAT64_MATH_RANDOM_TEST_TEMPLATE(TNAME, FNAME, LENGTH, LOW, HIGH) \
TEST(base_libm, TNAME) {                                                        \
    int len = (LENGTH);                                                         \
    std::vector<double> da(len);                                                \
    std::vector<double> ds1(len);                                               \
    std::vector<double> ds2(len);                                               \
    std::uniform_real_distribution<double> rnd(LOW, HIGH);                      \
    std::default_random_engine eng;                                             \
    for (int i = 0; i < len; ++i)                                               \
        da[i] = rnd(eng);                                                       \
    for (int i = 0; i < len; ++i) {                                             \
        ds1[i] = std::FNAME(da[i]);                                             \
        ds2[i] = ModuleBase::libm::FNAME(da[i]);                                \
    }                                                                           \
    for (int i = 0; i < len; i++)                                               \
        MY_EXPECT_DOUBLE_EQ(ds1[i], ds2[i]);                                       \
}                                                                               \

SINGLE_PARAM_FLOAT64_MATH_RANDOM_TEST_TEMPLATE(exp_random, exp, 50000, -1000, 1000)
SINGLE_PARAM_FLOAT64_MATH_RANDOM_TEST_TEMPLATE(sin_random, sin, 50000, -50000, 50000)
SINGLE_PARAM_FLOAT64_MATH_RANDOM_TEST_TEMPLATE(cos_random, cos, 50000, -50000, 50000)
SINGLE_PARAM_FLOAT64_MATH_RANDOM_TEST_TEMPLATE(sin_random2, sin, 100000, -DBL_MAX / M_PI, DBL_MAX / M_PI)
SINGLE_PARAM_FLOAT64_MATH_RANDOM_TEST_TEMPLATE(cos_random2, cos, 100000, -DBL_MAX / M_PI, DBL_MAX / M_PI)

#define SINGLE_PARAM_FLOAT64_MATH_SPEC_TEST_TEMPLATE(TNAME, FNAME, VALUE)       \
TEST(base_libm, TNAME) {                                                        \
    double ds1 = std::FNAME(VALUE);                                             \
    double ds2 = ModuleBase::libm::FNAME(VALUE);                                \
    MY_EXPECT_DOUBLE_EQ(ds1, ds2);                                              \
}                                                                               \

#define SINGLE_PARAM_FLOAT64_MATH_CORNER_TEST_TEMPLATE(FNAME)                   \
SINGLE_PARAM_FLOAT64_MATH_SPEC_TEST_TEMPLATE(FNAME ## _nz, FNAME, -0.0);        \
SINGLE_PARAM_FLOAT64_MATH_SPEC_TEST_TEMPLATE(FNAME ## _pz, FNAME, +0.0);        \
SINGLE_PARAM_FLOAT64_MATH_SPEC_TEST_TEMPLATE(FNAME ## _ninf, FNAME, -INFINITY); \
SINGLE_PARAM_FLOAT64_MATH_SPEC_TEST_TEMPLATE(FNAME ## _pinf, FNAME, +INFINITY); \
SINGLE_PARAM_FLOAT64_MATH_SPEC_TEST_TEMPLATE(FNAME ## _nnan, FNAME, -NAN);      \
SINGLE_PARAM_FLOAT64_MATH_SPEC_TEST_TEMPLATE(FNAME ## _pnan, FNAME, +NAN);      \
SINGLE_PARAM_FLOAT64_MATH_SPEC_TEST_TEMPLATE(FNAME ## _nmin, FNAME, -DBL_MIN);  \
SINGLE_PARAM_FLOAT64_MATH_SPEC_TEST_TEMPLATE(FNAME ## _nmax, FNAME, -DBL_MAX);  \
SINGLE_PARAM_FLOAT64_MATH_SPEC_TEST_TEMPLATE(FNAME ## _pmin, FNAME, +DBL_MIN);  \
SINGLE_PARAM_FLOAT64_MATH_SPEC_TEST_TEMPLATE(FNAME ## _pmax, FNAME, +DBL_MAX);  \
SINGLE_PARAM_FLOAT64_MATH_SPEC_TEST_TEMPLATE(FNAME ## _nmin2, FNAME, -DBL_MIN * DBL_MIN);\
SINGLE_PARAM_FLOAT64_MATH_SPEC_TEST_TEMPLATE(FNAME ## _nmax2, FNAME, -DBL_MAX * DBL_MAX);\
SINGLE_PARAM_FLOAT64_MATH_SPEC_TEST_TEMPLATE(FNAME ## _pmin2, FNAME, +DBL_MIN * DBL_MIN);\
SINGLE_PARAM_FLOAT64_MATH_SPEC_TEST_TEMPLATE(FNAME ## _pmax2, FNAME, +DBL_MAX * DBL_MAX);\
SINGLE_PARAM_FLOAT64_MATH_SPEC_TEST_TEMPLATE(FNAME ## _nlarge, FNAME, -DBL_MAX / M_PI);\
SINGLE_PARAM_FLOAT64_MATH_SPEC_TEST_TEMPLATE(FNAME ## _plarge, FNAME, +DBL_MAX / M_PI);\
SINGLE_PARAM_FLOAT64_MATH_SPEC_TEST_TEMPLATE(FNAME ## _nlarge2, FNAME, -105414350.0 * 2.0);\
SINGLE_PARAM_FLOAT64_MATH_SPEC_TEST_TEMPLATE(FNAME ## _plarge2, FNAME, +105414350.0 * 2.0);\

SINGLE_PARAM_FLOAT64_MATH_CORNER_TEST_TEMPLATE(sin)
SINGLE_PARAM_FLOAT64_MATH_CORNER_TEST_TEMPLATE(cos)
SINGLE_PARAM_FLOAT64_MATH_CORNER_TEST_TEMPLATE(exp)

TEST(base_libm, cexp_random)
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
        MY_EXPECT_COMPLEX_DOUBLE_EQ(ds1[i], ds2[i]);
    }
}

TEST(base_libm, cexp_spec)
{
    std::vector<std::complex<double>> da = {
        {+INFINITY, +0.0}, // 1
        {+INFINITY, -0.0},
        {-INFINITY, +0.0},
        {-INFINITY, -0.0},
        {+INFINITY, 2.0},  // 5
        {+INFINITY, 4.0},
        {-INFINITY, 2.0},
        {-INFINITY, 4.0},
        {+0.0, +INFINITY},
        {-0.0, +INFINITY}, // 10
        {+0.0, -INFINITY},
        {-0.0, -INFINITY},
        {+100.0, +INFINITY},
        {-100.0, +INFINITY},
        {+100.0, -INFINITY}, // 15
        {-100.0, -INFINITY},
        {+INFINITY, +INFINITY},
        {-INFINITY, +INFINITY},
        {+INFINITY, -INFINITY},
        {-INFINITY, -INFINITY}, // 20
        {+INFINITY, NAN},
        {-INFINITY, NAN},
        {+0.0, NAN},
        {+1.0, NAN},
        {NAN, +0.0}, // 25
        {NAN, -0.0},
        {NAN, +1.0},
        {NAN, +INFINITY},
        {NAN, NAN},
        {+0.0, -DBL_MIN * DBL_MIN}, // 30
        {+0.0, +DBL_MAX * DBL_MAX},
        {+0.0, -DBL_MIN * DBL_MIN},
        {+0.0, +DBL_MAX * DBL_MAX},
        {+INFINITY, -DBL_MIN * DBL_MIN},
        {+INFINITY, +DBL_MAX * DBL_MAX}, // 35
        {+INFINITY, -DBL_MIN * DBL_MIN},
        {+INFINITY, +DBL_MAX * DBL_MAX},
        {+0.0, -DBL_MIN},
        {+0.0, +DBL_MAX},
        {+0.0, -DBL_MIN}, // 40
        {+0.0, +DBL_MAX},
        {+INFINITY, -DBL_MIN},
        {+INFINITY, +DBL_MAX},
        {+INFINITY, -DBL_MIN},
        {+INFINITY, +DBL_MAX}, // 45
        {+DBL_MAX, +0.0},
        {-DBL_MAX, +0.0},
        {+DBL_MAX / M_PI, +0.0},
        {-DBL_MAX / M_PI, +0.0},
        {+DBL_MIN, +0.0}, // 50
        {-DBL_MIN, +0.0},
    };

    int len = da.size();
    std::vector<std::complex<double>> ds1(len);
    std::vector<std::complex<double>> ds2(len);

    for (int i = 0; i < len; ++i) {
        ds1[i] = std::exp(da[i]);
        ds2[i] = ModuleBase::libm::exp(da[i]);
    }

    for (int i = 0; i < len; i++) {
        MY_EXPECT_COMPLEX_DOUBLE_EQ(ds1[i], ds2[i]);
    }
}
