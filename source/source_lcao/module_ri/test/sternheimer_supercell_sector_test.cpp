#include "source_lcao/module_ri/sternheimer_supercell_sector.h"

#include <cmath>
#include <complex>
#include <gtest/gtest.h>
#include <vector>

namespace
{

using Complex = std::complex<double>;

std::vector<Complex> sector_vector(const int primitive_orbital, const double k)
{
    const double inverse_sqrt_two = 1.0 / std::sqrt(2.0);
    const Complex phase = std::exp(Complex(0.0, 2.0 * std::acos(-1.0) * k));
    std::vector<Complex> result(4, Complex(0.0, 0.0));
    result[static_cast<std::size_t>(primitive_orbital)] = inverse_sqrt_two;
    result[static_cast<std::size_t>(2 + primitive_orbital)] = phase * inverse_sqrt_two;
    return result;
}

std::vector<Complex> mixed_sector_vector(const std::array<Complex, 2>& primitive_coefficients,
                                         const double k)
{
    const double inverse_sqrt_two = 1.0 / std::sqrt(2.0);
    const Complex phase = std::exp(Complex(0.0, 2.0 * std::acos(-1.0) * k));
    return {inverse_sqrt_two * primitive_coefficients[0],
            inverse_sqrt_two * primitive_coefficients[1],
            phase * inverse_sqrt_two * primitive_coefficients[0],
            phase * inverse_sqrt_two * primitive_coefficients[1]};
}

} // namespace

TEST(SternheimerSupercellSector, RecoversSelectedTranslationSectorFromCompleteEigensystem)
{
    const std::vector<double> eigenvalues{1.0, 2.0, 4.0, 7.0};
    const std::vector<std::vector<Complex>> eigenvectors{
        sector_vector(0, 0.0),
        sector_vector(1, 0.0),
        sector_vector(0, 0.5),
        sector_vector(1, 0.5),
    };

    const auto sector = ModuleRI::recover_sternheimer_supercell_sector(
        eigenvalues, eigenvectors, {2, 1, 1}, {0.5, 0.0, 0.0});

    ASSERT_EQ(sector.eigenvalues.size(), 2U);
    EXPECT_NEAR(sector.eigenvalues[0], 4.0, 1.0e-12);
    EXPECT_NEAR(sector.eigenvalues[1], 7.0, 1.0e-12);
    ASSERT_EQ(sector.coefficients.size(), 2U);
    EXPECT_LT(sector.max_orthonormality_error, 1.0e-12);
    EXPECT_LT(sector.max_full_space_residual, 1.0e-12);
    for (const auto& coefficients: sector.coefficients)
    {
        ASSERT_EQ(coefficients.size(), 4U);
        EXPECT_NEAR(coefficients[0].real(), -coefficients[2].real(), 1.0e-12);
        EXPECT_NEAR(coefficients[1].real(), -coefficients[3].real(), 1.0e-12);
    }
}

TEST(SternheimerSupercellSector, RejectsIncompleteSquareEigensystem)
{
    const std::vector<double> eigenvalues{1.0, 2.0, 3.0};
    const std::vector<std::vector<Complex>> eigenvectors{
        std::vector<Complex>(4, Complex(0.0, 0.0)),
        std::vector<Complex>(4, Complex(0.0, 0.0)),
        std::vector<Complex>(4, Complex(0.0, 0.0)),
    };
    EXPECT_THROW(ModuleRI::recover_sternheimer_supercell_sector(
                     eigenvalues, eigenvectors, {2, 1, 1}, {0.0, 0.0, 0.0}),
                 std::invalid_argument);
}

TEST(SternheimerSupercellSector, RecoversDenseComplexNonunitarySectorWithoutLayoutAmbiguity)
{
    const std::vector<double> eigenvalues{1.0, 3.0, 5.0, 8.0};
    const std::vector<std::vector<Complex>> eigenvectors{
        mixed_sector_vector({Complex(1.0, 0.25), Complex(0.4, -0.3)}, 0.0),
        mixed_sector_vector({Complex(-0.2, 0.7), Complex(1.3, 0.1)}, 0.0),
        mixed_sector_vector({Complex(0.8, -0.4), Complex(0.1, 0.9)}, 0.5),
        mixed_sector_vector({Complex(-0.6, 0.2), Complex(1.1, -0.5)}, 0.5),
    };

    const auto sector = ModuleRI::recover_sternheimer_supercell_sector(
        eigenvalues, eigenvectors, {2, 1, 1}, {0.5, 0.0, 0.0});

    ASSERT_EQ(sector.eigenvalues.size(), 2U);
    EXPECT_NEAR(sector.eigenvalues[0], 5.0, 1.0e-11);
    EXPECT_NEAR(sector.eigenvalues[1], 8.0, 1.0e-11);
    EXPECT_LT(sector.max_orthonormality_error, 1.0e-11);
    EXPECT_LT(sector.max_full_space_residual, 1.0e-11);
    for (const auto& coefficients: sector.coefficients)
    {
        ASSERT_EQ(coefficients.size(), 4U);
        EXPECT_NEAR(std::abs(coefficients[0] + coefficients[2]), 0.0, 1.0e-11);
        EXPECT_NEAR(std::abs(coefficients[1] + coefficients[3]), 0.0, 1.0e-11);
    }
}
