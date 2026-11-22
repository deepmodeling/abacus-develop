#include <limits>
#include <stdexcept>

#include <gtest/gtest.h>

#include "source_base/matrix.h"

// Issue #8342: NLOCAL = 377884 makes NLOCAL^2 about 1.4e11, far above
// INT_MAX.  The helper must keep the product in size_t.
TEST(MatrixSize, LargeLcaoDimensionsDoNotOverflow)
{
    const std::size_t size = ModuleBase::matrix_size(377884, 377884);
    const std::size_t expected = static_cast<std::size_t>(377884) * 377884;
    EXPECT_EQ(size, expected);
    EXPECT_GT(size, static_cast<std::size_t>(std::numeric_limits<int>::max()));
}

TEST(MatrixSize, ZeroDimension)
{
    EXPECT_EQ(ModuleBase::matrix_size(0, 123), 0u);
    EXPECT_EQ(ModuleBase::matrix_size(123, 0), 0u);
}

TEST(MatrixSize, NegativeDimensionThrows)
{
    EXPECT_THROW(ModuleBase::matrix_size(-1, 10), std::length_error);
    EXPECT_THROW(ModuleBase::matrix_size(10, -1), std::length_error);
}

// A row-major index can overflow int even when the total element count is
// still small enough for some allocators, so indexing is separate from sizing.
TEST(MatrixIndex, LargeRowIndexDoesNotOverflow)
{
    const std::size_t index = ModuleBase::matrix_index(100000, 200000, 377884);
    const std::size_t expected = static_cast<std::size_t>(100000) * 377884 + 200000;
    EXPECT_EQ(index, expected);
}

TEST(MatrixIndex, RowMajorOrder)
{
    EXPECT_EQ(ModuleBase::matrix_index(0, 0, 4), 0u);
    EXPECT_EQ(ModuleBase::matrix_index(0, 3, 4), 3u);
    EXPECT_EQ(ModuleBase::matrix_index(2, 1, 4), 9u);
}

TEST(Matrix, ElementAccessUsesSizeTIndex)
{
    ModuleBase::matrix m(3, 4, true);
    m(2, 3) = 5.0;
    EXPECT_DOUBLE_EQ(m(2, 3), 5.0);
}

TEST(Matrix, NegativeConstructorThrows)
{
    EXPECT_THROW(ModuleBase::matrix m(-1, 4), std::length_error);
}

TEST(Matrix, CreateNegativeThrows)
{
    ModuleBase::matrix m;
    EXPECT_THROW(m.create(-1, 4), std::length_error);
}

TEST(Matrix, EmptyMatrixHasNoStorage)
{
    ModuleBase::matrix m(0, 4);
    EXPECT_EQ(m.nr, 0);
    EXPECT_EQ(m.nc, 4);
    EXPECT_EQ(m.c, nullptr);
}
