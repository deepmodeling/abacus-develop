#include "../rpa_lri_detail.h"

#include <RI/global/Tensor.h>
#include <gtest/gtest.h>

#include <complex>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace
{
using RpaLriDetail::checked_add_u64;
using RpaLriDetail::checked_i32_from_int;
using RpaLriDetail::checked_i32_from_size;
using RpaLriDetail::checked_i64_from_u64;
using RpaLriDetail::checked_mul_u64;
using RpaLriDetail::coulomb_atom_pair_index;
using RpaLriDetail::has_valid_matrix_shape;
using RpaLriDetail::real_as_double;
using RpaLriDetail::sum_int_vector;
}

TEST(RpaLriDetail, CoulombAtomPairIndexUsesUpperTriangleOrder)
{
    const std::size_t natoms = 4;
    const std::vector<std::size_t> expected = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    std::size_t index = 0;
    for (std::size_t I = 0; I < natoms; ++I)
    {
        for (std::size_t J = I; J < natoms; ++J)
        {
            EXPECT_EQ(coulomb_atom_pair_index(I, J, natoms), expected[index]);
            ++index;
        }
    }
    EXPECT_THROW(coulomb_atom_pair_index(2, 1, natoms), std::runtime_error);
}

TEST(RpaLriDetail, CheckedUnsignedArithmeticRejectsOverflow)
{
    const unsigned long long max_value = std::numeric_limits<unsigned long long>::max();
    EXPECT_EQ(checked_mul_u64(6, 7, "product"), 42);
    EXPECT_EQ(checked_add_u64(6, 7, "sum"), 13);
    EXPECT_THROW(checked_mul_u64(max_value, 2, "product"), std::runtime_error);
    EXPECT_THROW(checked_add_u64(max_value, 1, "sum"), std::runtime_error);
}

TEST(RpaLriDetail, CheckedIntegerConversionsRejectOverflowAndNegativeValues)
{
    EXPECT_EQ(checked_i64_from_u64(17, "value"), 17);
    EXPECT_EQ(checked_i32_from_size(17, "value"), 17);
    EXPECT_EQ(checked_i32_from_int(17, "value"), 17);
    EXPECT_THROW(checked_i64_from_u64(std::numeric_limits<unsigned long long>::max(), "value"),
                 std::runtime_error);
    EXPECT_THROW(checked_i32_from_size(std::numeric_limits<std::size_t>::max(), "value"), std::runtime_error);
    EXPECT_THROW(checked_i32_from_int(-1, "value"), std::runtime_error);
}

TEST(RpaLriDetail, SumIntVectorRejectsPositiveOverflow)
{
    EXPECT_EQ(sum_int_vector({1, 2, 3}), 6);
    EXPECT_EQ(sum_int_vector({}), 0);
    EXPECT_THROW(sum_int_vector({std::numeric_limits<int>::max(), 1}), std::runtime_error);
}

TEST(RpaLriDetail, RealPartConversionAndMatrixShapeValidation)
{
    EXPECT_DOUBLE_EQ(real_as_double(2.5), 2.5);
    EXPECT_DOUBLE_EQ(real_as_double(std::complex<double>(2.5, -4.0)), 2.5);

    RI::Tensor<double> tensor;
    EXPECT_FALSE(has_valid_matrix_shape(tensor));
    tensor.shape = {2, 3};
    EXPECT_TRUE(has_valid_matrix_shape(tensor));
    tensor.shape = {2};
    EXPECT_FALSE(has_valid_matrix_shape(tensor));
    tensor.shape = {0, 3};
    EXPECT_FALSE(has_valid_matrix_shape(tensor));
}
