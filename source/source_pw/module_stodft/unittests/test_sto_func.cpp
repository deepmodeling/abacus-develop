#include "source_estate/occupy.h"
#include "source_pw/module_stodft/sto_func.h"

#include <cmath>
#include <gtest/gtest.h>
#include <limits>
#include <type_traits>

// Supply only the external smearing setting consumed by the real constructor.
// All Sto_Func methods are linked from sto_func.cpp, without the SCF machinery.
double Occupy::gaussian_parameter = 0.125;

namespace
{
template <typename Real> void expect_relative(Real actual, long double expected)
{
    const long double magnitude = std::abs(expected);
    const long double tolerance = 32 * std::numeric_limits<Real>::epsilon() * magnitude;
    const long double observed = actual;
    // A relative tolerance is essential: returning zero must fail even for tails.
    EXPECT_NEAR(observed, expected, tolerance);
}

long double fermi_reference(long double x)
{
    const long double negative_x = -x;
    const long double boltzmann = std::exp(negative_x);
    return boltzmann / (1 + boltzmann);
}

template <typename Real> class StoFuncTest : public testing::Test
{
  protected:
    Sto_Func<Real> function;
    Real lower = -3;
    Real upper = 5;

    void SetUp() override
    {
        function.set_E_range(&lower, &upper);
        function.tem = 2;
        function.mu = 1;
        function.sigma = 2;
        function.targ_e = 1;
        function.t = 0;
    }

    void use_identity_range()
    {
        lower = -1;
        upper = 1;
        function.tem = 1;
        function.mu = 0;
        function.sigma = 1;
        function.targ_e = 0;
    }
};

#ifdef __FLOAT_FFTW
using Precisions = testing::Types<double, float>;
#else
using Precisions = testing::Types<double>;
#endif
TYPED_TEST_SUITE(StoFuncTest, Precisions);

TYPED_TEST(StoFuncTest, ConstructorUsesSmearingAndInitializesControls)
{
    const Sto_Func<TypeParam> fresh;
    const TypeParam smearing = Occupy::gaussian_parameter;
    EXPECT_EQ(fresh.tem, smearing);
    EXPECT_EQ(fresh.mu, 0);
    EXPECT_EQ(fresh.t, 0);
    EXPECT_EQ(fresh.sigma, 0);
    EXPECT_EQ(fresh.targ_e, 0);
    EXPECT_EQ(fresh.Emin, nullptr);
    EXPECT_EQ(fresh.Emax, nullptr);
}

TYPED_TEST(StoFuncTest, FermiValuesAndRootsMatchIndependentReference)
{
    const TypeParam scaled[] = {-1,
        -0.5,
        0,
        0.5,
        1};
    const TypeParam energies[] = {-3,
        -1,
        1,
        3,
        5};
    for (int i = 0; i < 5; ++i)
    {
        SCOPED_TRACE(i);
        const long double x = i - 2;
        const long double negative_x = -x;
        const long double occupation = fermi_reference(x);
        const long double hole = fermi_reference(negative_x);
        const long double root = std::sqrt(occupation);
        const long double root_hole = std::sqrt(hole);
        const long double weighted = energies[i] * occupation;
        const TypeParam direct = this->function.fd(energies[i]);
        const TypeParam normalized = this->function.nfd(scaled[i]);
        const TypeParam direct_root = this->function.root_fd(energies[i]);
        const TypeParam normalized_root = this->function.nroot_fd(scaled[i]);
        const TypeParam normalized_hole = this->function.nroot_mfd(scaled[i]);
        const TypeParam energy_weighted = this->function.nxfd(scaled[i]);
        expect_relative(direct, occupation);
        expect_relative(normalized, occupation);
        expect_relative(direct_root, root);
        expect_relative(normalized_root, root);
        expect_relative(normalized_hole, root_hole);
        expect_relative(energy_weighted, weighted);
    }
}

TYPED_TEST(StoFuncTest, FermiCutoffKeepsTheBoundaryAndRejectsItsNextNeighbor)
{
    this->use_identity_range();
    const TypeParam boundary = 36;
    const TypeParam infinity = std::numeric_limits<TypeParam>::infinity();
    const TypeParam below = std::nextafter(boundary, TypeParam(0));
    const TypeParam above = std::nextafter(boundary, infinity);
    const TypeParam inside[] = {below, boundary};
    for (const TypeParam energy : inside)
    {
        SCOPED_TRACE(energy);
        const long double expected = fermi_reference(energy);
        const long double weighted = energy * expected;
        const TypeParam direct = this->function.fd(energy);
        const TypeParam normalized = this->function.nfd(energy);
        const TypeParam energy_weighted = this->function.nxfd(energy);
        expect_relative(direct, expected);
        expect_relative(normalized, expected);
        expect_relative(energy_weighted, weighted);
    }
    const TypeParam direct = this->function.fd(above);
    const TypeParam normalized = this->function.nfd(above);
    const TypeParam weighted = this->function.nxfd(above);
    EXPECT_EQ(direct, 0);
    EXPECT_EQ(normalized, 0);
    EXPECT_EQ(weighted, 0);
}

TYPED_TEST(StoFuncTest, RootCutoffsKeepTheirBoundariesOnOppositeTails)
{
    this->use_identity_range();
    const TypeParam boundary = 72;
    const TypeParam infinity = std::numeric_limits<TypeParam>::infinity();
    const TypeParam below = std::nextafter(boundary, TypeParam(0));
    const TypeParam above = std::nextafter(boundary, infinity);
    const TypeParam inside[] = {below, boundary};
    for (const TypeParam energy : inside)
    {
        SCOPED_TRACE(energy);
        const long double occupation = fermi_reference(energy);
        const long double expected = std::sqrt(occupation);
        const TypeParam negative_energy = -energy;
        const TypeParam direct = this->function.root_fd(energy);
        const TypeParam normalized = this->function.nroot_fd(energy);
        const TypeParam hole = this->function.nroot_mfd(negative_energy);
        expect_relative(direct, expected);
        expect_relative(normalized, expected);
        expect_relative(hole, expected);
    }
    const TypeParam negative_above = -above;
    const TypeParam direct = this->function.root_fd(above);
    const TypeParam normalized = this->function.nroot_fd(above);
    const TypeParam hole = this->function.nroot_mfd(negative_above);
    EXPECT_EQ(direct, 0);
    EXPECT_EQ(normalized, 0);
    EXPECT_EQ(hole, 0);

    const TypeParam intermediate = 40;
    const long double occupation = fermi_reference(intermediate);
    const long double expected = std::sqrt(occupation);
    const TypeParam root = this->function.nroot_fd(intermediate);
    const TypeParam truncated = this->function.nfd(intermediate);
    expect_relative(root, expected);
    EXPECT_EQ(truncated, 0);
}

TYPED_TEST(StoFuncTest, EntropyMatchesHyperbolicReferenceAndHasCorrectSign)
{
    const TypeParam scaled[] = {-1,
        -0.5,
        0,
        0.5,
        1};
    const TypeParam energies[] = {-3,
        -1,
        1,
        3,
        5};
    for (int i = 0; i < 5; ++i)
    {
        SCOPED_TRACE(i);
        const long double half_x = (i - 2) / 2.0L;
        const long double twice_cosh = 2 * std::cosh(half_x);
        const long double expected = half_x * std::tanh(half_x) - std::log(twice_cosh);
        const long double positive_entropy = -expected;
        const long double root_expected = std::sqrt(positive_entropy);
        const TypeParam direct = this->function.fdlnfd(energies[i]);
        const TypeParam normalized = this->function.nfdlnfd(scaled[i]);
        const TypeParam root = this->function.n_root_fdlnfd(scaled[i]);
        expect_relative(direct, expected);
        expect_relative(normalized, expected);
        expect_relative(root, root_expected);
        EXPECT_LT(normalized, 0);
        EXPECT_GT(root, 0);
    }
}

TYPED_TEST(StoFuncTest, EntropyTailGuardsAvoidLogOfZero)
{
    this->use_identity_range();
    const TypeParam boundary = 36;
    const TypeParam infinity = std::numeric_limits<TypeParam>::infinity();
    const TypeParam above = std::nextafter(boundary, infinity);
    const TypeParam negative_above = -above;
    const TypeParam outside[] = {negative_above, above};
    for (const TypeParam energy : outside)
    {
        SCOPED_TRACE(energy);
        const TypeParam direct = this->function.fdlnfd(energy);
        const TypeParam normalized = this->function.nfdlnfd(energy);
        const TypeParam root = this->function.n_root_fdlnfd(energy);
        EXPECT_EQ(direct, 0);
        EXPECT_EQ(normalized, 0);
        EXPECT_EQ(root, 0);
    }
    const TypeParam direct_boundary = this->function.fdlnfd(boundary);
    const TypeParam normalized_boundary = this->function.nfdlnfd(boundary);
    const TypeParam root_boundary = this->function.n_root_fdlnfd(boundary);
    EXPECT_LT(direct_boundary, 0);
    EXPECT_LT(normalized_boundary, 0);
    EXPECT_GT(root_boundary, 0);

    // Inside the cutoff, float rounds the occupation to one; double does not.
    const TypeParam negative_boundary = -boundary;
    const TypeParam direct = this->function.fdlnfd(negative_boundary);
    const TypeParam normalized = this->function.nfdlnfd(negative_boundary);
    const TypeParam root = this->function.n_root_fdlnfd(negative_boundary);
    if (std::is_same<TypeParam, float>::value)
    {
        EXPECT_EQ(direct, 0);
        EXPECT_EQ(normalized, 0);
        EXPECT_EQ(root, 0);
    }
    else
    {
        EXPECT_LT(direct, 0);
        EXPECT_LT(normalized, 0);
        EXPECT_GT(root, 0);
    }
}

TYPED_TEST(StoFuncTest, TrigonometricFunctionsUseMappedEnergyAndTime)
{
    const TypeParam scaled[] = {-1,
        -0.5,
        0,
        0.5,
        1};
    const TypeParam energies[] = {-3,
        -1,
        1,
        3,
        5};
    this->function.t = 0.25;
    for (int i = 0; i < 5; ++i)
    {
        SCOPED_TRACE(i);
        const long double phase = energies[i] * 0.25L;
        const long double cosine = std::cos(phase);
        const long double sine = std::sin(phase);
        const long double negative_sine = -sine;
        const TypeParam actual_cosine = this->function.ncos(scaled[i]);
        const TypeParam actual_sine = this->function.nsin(scaled[i]);
        const TypeParam actual_negative = this->function.n_sin(scaled[i]);
        expect_relative(actual_cosine, cosine);
        expect_relative(actual_sine, sine);
        expect_relative(actual_negative, negative_sine);
    }
    this->function.t = 0;
    const TypeParam cosine = this->function.ncos(1);
    const TypeParam sine = this->function.nsin(1);
    const TypeParam negative_sine = this->function.n_sin(1);
    EXPECT_EQ(cosine, 1);
    EXPECT_EQ(sine, 0);
    EXPECT_EQ(negative_sine, 0);
}

TYPED_TEST(StoFuncTest, GaussiansMatchNormalDistributionAndWidthScaling)
{
    const long double pi = std::acos(-1.0L);
    const long double two_pi = 2 * pi;
    const long double normalization = std::sqrt(two_pi);
    const TypeParam widths[] = {0.5, 2};
    const TypeParam offsets[] = {-2,
        -1,
        0,
        1,
        2};
    for (const TypeParam width : widths)
    {
        SCOPED_TRACE(width);
        this->function.sigma = width;
        for (const TypeParam offset : offsets)
        {
            SCOPED_TRACE(offset);
            const TypeParam energy = 1 + width * offset;
            const TypeParam scaled = (energy - 1) / 4;
            const long double exponent = -0.5L * offset * offset;
            const long double expected = std::exp(exponent) / (normalization * width);
            const long double root_expected = std::sqrt(expected);
            const TypeParam direct = this->function.gauss(energy);
            const TypeParam normalized = this->function.ngauss(scaled);
            const TypeParam root = this->function.nroot_gauss(scaled);
            expect_relative(direct, expected);
            expect_relative(normalized, expected);
            expect_relative(root, root_expected);
        }
    }
}

TYPED_TEST(StoFuncTest, GaussianCutoffKeepsExponentSeventyTwo)
{
    this->use_identity_range();
    const long double pi = std::acos(-1.0L);
    const long double two_pi = 2 * pi;
    const long double normalization = std::sqrt(two_pi);
    const TypeParam boundary = 12; // E^2 / (2 sigma^2) == 72 exactly.
    const TypeParam infinity = std::numeric_limits<TypeParam>::infinity();
    const TypeParam below = std::nextafter(boundary, TypeParam(0));
    const TypeParam above = std::nextafter(boundary, infinity);
    const TypeParam signs[] = {-1, 1};
    for (const TypeParam sign : signs)
    {
        const TypeParam inside[] = {sign * below, sign * boundary};
        for (const TypeParam energy : inside)
        {
            SCOPED_TRACE(energy);
            const long double exponent = -0.5L * energy * energy;
            const long double expected = std::exp(exponent) / normalization;
            const TypeParam direct = this->function.gauss(energy);
            const TypeParam normalized = this->function.ngauss(energy);
            // Rounding the stored REAL exponent is amplified by |exponent|
            // in exp(). Account for that conditioning near a = 72, while
            // keeping a relative tolerance that still rejects a zero tail.
            const long double exponent_magnitude = std::abs(exponent);
            const long double relative_error = (2 + exponent_magnitude) * std::numeric_limits<TypeParam>::epsilon();
            const long double tolerance = relative_error * expected;
            const long double direct_wide = direct;
            const long double normalized_wide = normalized;
            EXPECT_NEAR(direct_wide, expected, tolerance);
            EXPECT_NEAR(normalized_wide, expected, tolerance);
        }
        const TypeParam outside = sign * above;
        const TypeParam direct = this->function.gauss(outside);
        const TypeParam normalized = this->function.ngauss(outside);
        EXPECT_EQ(direct, 0);
        EXPECT_EQ(normalized, 0);
    }
}

TYPED_TEST(StoFuncTest, RootGaussianHasItsOwnWiderTailCutoff)
{
    this->use_identity_range();
    const long double pi = std::acos(-1.0L);
    const long double two_pi = 2 * pi;
    const long double normalization = std::sqrt(two_pi);
    const TypeParam inside[] = {-16,
        -14,
        14,
        16};
    for (const TypeParam energy : inside)
    {
        SCOPED_TRACE(energy);
        const long double exponent = -0.5L * energy * energy;
        const long double gaussian = std::exp(exponent) / normalization;
        const long double expected = std::sqrt(gaussian);
        const TypeParam root = this->function.nroot_gauss(energy);
        const TypeParam truncated = this->function.ngauss(energy);
        expect_relative(root, expected);
        EXPECT_EQ(truncated, 0);
    }
    const TypeParam outside[] = {-18, 18};
    for (const TypeParam energy : outside)
    {
        const TypeParam root = this->function.nroot_gauss(energy);
        EXPECT_EQ(root, 0);
    }
}

TYPED_TEST(StoFuncTest, EnergyRangeCanBeReboundAndUpdatedInPlace)
{
    TypeParam new_lower = -7;
    TypeParam new_upper = 1;
    this->function.set_E_range(&new_lower, &new_upper);
    // Mutating the old range must not affect the new binding.
    this->lower = 20;
    this->upper = 30;
    const TypeParam first = this->function.nfd(0);
    const long double expected_first = fermi_reference(-2);
    expect_relative(first, expected_first);

    new_lower = -3;
    new_upper = 5;
    const TypeParam second = this->function.nfd(0);
    EXPECT_EQ(second, TypeParam(0.5));
    const TypeParam last_energy = this->function.nxfd(1);
    const long double occupation = fermi_reference(2);
    const long double expected_last = 5 * occupation;
    expect_relative(last_energy, expected_last);
}
} // namespace
