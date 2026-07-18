#include "source_lcao/module_ri/sternheimer_siab_data.h"
#include "source_lcao/module_ri/sternheimer_siab_overlap.h"

#include <complex>
#include <gtest/gtest.h>
#include <limits>
#include <stdexcept>
#include <vector>

namespace
{

using Complex = std::complex<double>;
using GridVector = std::vector<Complex>;
using PrimitiveGrid = std::vector<GridVector>;
namespace siab = module_ri::sternheimer_siab;

const GridVector reference_wavefunction = {
    {1.0, 1.0},
    {2.0, 0.0},
    {0.0, -1.0},
    {1.0, 0.0},
};

const PrimitiveGrid primitives = {
    {{1.0, 0.0}, {0.0, 0.0}, {1.0, 0.0}, {0.0, 0.0}},
    {{0.0, 0.0}, {1.0, 0.0}, {0.0, 0.0}, {1.0, 0.0}},
};

} // namespace

TEST(SternheimerSIABData, ExposesVersionOneFields)
{
    const siab::PrimitiveBlock block{"H", 0, 1, -1, 3, 6};
    EXPECT_EQ(block.element, "H");
    EXPECT_EQ(block.atom_index, 0);
    EXPECT_EQ(block.l, 1);
    EXPECT_EQ(block.m, -1);
    EXPECT_EQ(block.n_primitive, 3);
    EXPECT_EQ(block.offset, 6);

    const siab::ReferenceRow row{2, 4, 0.5, 1.0, 0.25, 1.5, {{0.75, -0.5}}};
    EXPECT_EQ(row.occupied_state, 2);
    EXPECT_EQ(row.auxiliary_channel, 4);
    EXPECT_DOUBLE_EQ(row.frequency_ha, 0.5);
    EXPECT_DOUBLE_EQ(row.occupation, 1.0);
    EXPECT_DOUBLE_EQ(row.frequency_weight, 0.25);
    EXPECT_DOUBLE_EQ(row.norm, 1.5);
    ASSERT_EQ(row.q.size(), 1);
    EXPECT_EQ(row.q[0], Complex(0.75, -0.5));
}

TEST(SternheimerSIABOverlap, MatchesPlanFourPointFixture)
{
    constexpr double delta_omega = 0.25;

    EXPECT_NEAR(siab::norm(reference_wavefunction, delta_omega), 2.0, 1.0e-14);

    const auto q = siab::overlap_q(reference_wavefunction, primitives, delta_omega);
    ASSERT_EQ(q.size(), 2);
    EXPECT_NEAR(q[0].real(), 0.25, 1.0e-14);
    EXPECT_NEAR(q[0].imag(), 0.00, 1.0e-14);
    EXPECT_NEAR(q[1].real(), 0.75, 1.0e-14);
    EXPECT_NEAR(q[1].imag(), 0.00, 1.0e-14);
}

TEST(SternheimerSIABOverlap, ConjugatesReferenceWavefunctionInQ)
{
    const GridVector complex_reference = {{1.0, 1.0}, {0.0, 0.0}};
    const PrimitiveGrid real_primitive = {{{1.0, 0.0}, {0.0, 0.0}}};

    const auto q = siab::overlap_q(complex_reference, real_primitive, 0.5);

    ASSERT_EQ(q.size(), 1);
    EXPECT_NEAR(q[0].real(), 0.5, 1.0e-14);
    EXPECT_NEAR(q[0].imag(), -0.5, 1.0e-14);
}

TEST(SternheimerSIABOverlap, FlattensPrimitiveOverlapInRowMajorOrder)
{
    const auto s = siab::overlap_s(primitives, 0.25);

    ASSERT_EQ(s.size(), 4);
    EXPECT_NEAR(s[0].real(), 0.5, 1.0e-14);
    EXPECT_NEAR(s[0].imag(), 0.0, 1.0e-14);
    EXPECT_NEAR(s[1].real(), 0.0, 1.0e-14);
    EXPECT_NEAR(s[1].imag(), 0.0, 1.0e-14);
    EXPECT_NEAR(s[2].real(), 0.0, 1.0e-14);
    EXPECT_NEAR(s[2].imag(), 0.0, 1.0e-14);
    EXPECT_NEAR(s[3].real(), 0.5, 1.0e-14);
    EXPECT_NEAR(s[3].imag(), 0.0, 1.0e-14);

    constexpr std::size_t nprimitive = 2;
    for (std::size_t i = 0; i != nprimitive; ++i)
    {
        for (std::size_t j = 0; j != nprimitive; ++j)
        {
            EXPECT_NEAR(s[i * nprimitive + j].real(), s[j * nprimitive + i].real(), 1.0e-14);
            EXPECT_NEAR(s[i * nprimitive + j].imag(), -s[j * nprimitive + i].imag(), 1.0e-14);
        }
    }
}

TEST(SternheimerSIABOverlap, ProducesComplexHermitianPrimitiveOverlap)
{
    const PrimitiveGrid complex_primitives = {
        {{1.0, 1.0}, {0.0, 0.0}},
        {{1.0, 0.0}, {0.0, 0.0}},
    };

    const auto s = siab::overlap_s(complex_primitives, 0.5);

    ASSERT_EQ(s.size(), 4);
    EXPECT_NEAR(s[0].real(), 1.0, 1.0e-14);
    EXPECT_NEAR(s[0].imag(), 0.0, 1.0e-14);
    EXPECT_NEAR(s[1].real(), 0.5, 1.0e-14);
    EXPECT_NEAR(s[1].imag(), -0.5, 1.0e-14);
    EXPECT_NEAR(s[2].real(), 0.5, 1.0e-14);
    EXPECT_NEAR(s[2].imag(), 0.5, 1.0e-14);
    EXPECT_NEAR(s[3].real(), 0.5, 1.0e-14);
    EXPECT_NEAR(s[3].imag(), 0.0, 1.0e-14);
}

TEST(SternheimerSIABOverlap, RejectsInvalidGridVolume)
{
    EXPECT_THROW(siab::norm(reference_wavefunction, 0.0), std::invalid_argument);
    EXPECT_THROW(siab::overlap_q(reference_wavefunction, primitives, -0.25), std::invalid_argument);
    EXPECT_THROW(siab::overlap_s(primitives, std::numeric_limits<double>::infinity()), std::invalid_argument);
    EXPECT_THROW(siab::norm(reference_wavefunction, std::numeric_limits<double>::quiet_NaN()),
                 std::invalid_argument);
}

TEST(SternheimerSIABOverlap, EmptyLocalSlabContributesZero)
{
    const GridVector empty_grid;
    const PrimitiveGrid empty_primitives(2);

    EXPECT_DOUBLE_EQ(siab::norm(empty_grid, 0.25), 0.0);

    const auto q = siab::overlap_q(empty_grid, empty_primitives, 0.25);
    ASSERT_EQ(q.size(), 2);
    EXPECT_EQ(q[0], Complex(0.0, 0.0));
    EXPECT_EQ(q[1], Complex(0.0, 0.0));

    const auto s = siab::overlap_s(empty_primitives, 0.25);
    ASSERT_EQ(s.size(), 4);
    for (const Complex& value : s)
    {
        EXPECT_EQ(value, Complex(0.0, 0.0));
    }
}

TEST(SternheimerSIABOverlap, RejectsMissingBasisAndMismatchedGrids)
{
    EXPECT_THROW(siab::overlap_q(reference_wavefunction, {}, 0.25), std::invalid_argument);
    EXPECT_THROW(siab::overlap_s({}, 0.25), std::invalid_argument);

    PrimitiveGrid mismatched = primitives;
    mismatched[1].pop_back();
    EXPECT_THROW(siab::overlap_q(reference_wavefunction, mismatched, 0.25), std::invalid_argument);
    EXPECT_THROW(siab::overlap_s(mismatched, 0.25), std::invalid_argument);

    const PrimitiveGrid wrong_reference_size = {{{1.0, 0.0}, {0.0, 0.0}}};
    EXPECT_THROW(siab::overlap_q(reference_wavefunction, wrong_reference_size, 0.25), std::invalid_argument);

    const PrimitiveGrid mixed_empty_and_nonempty = {{}, {{1.0, 0.0}}};
    EXPECT_THROW(siab::overlap_q({}, mixed_empty_and_nonempty, 0.25), std::invalid_argument);
    EXPECT_THROW(siab::overlap_s(mixed_empty_and_nonempty, 0.25), std::invalid_argument);
}
