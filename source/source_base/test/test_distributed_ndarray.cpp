#ifdef NDEBUG
#undef NDEBUG
#endif

#include "source_base/module_container/NDArray/distributed_ndarray.h"

#include "gtest/gtest.h"

#include <cstddef>
#include <memory>
#include <type_traits>

static_assert(!std::is_copy_constructible<DistributedNDArray<double>>::value,
              "DistributedNDArray must not be copy constructible");
static_assert(!std::is_copy_assignable<DistributedNDArray<double>>::value,
              "DistributedNDArray must not be copy assignable");
static_assert(std::is_nothrow_move_constructible<DistributedNDArray<double>>::value,
              "DistributedNDArray move construction must be noexcept");
static_assert(std::is_nothrow_move_assignable<DistributedNDArray<double>>::value,
              "DistributedNDArray move assignment must be noexcept");

namespace
{
std::shared_ptr<BlockCyclicDistributionInformation> make_serial_distribution(
    const int global_rows,
    const int global_cols)
{
    std::shared_ptr<BlockCyclicDistributionInformation> distribution
        = std::make_shared<BlockCyclicDistributionInformation>();
    distribution->set_serial(global_rows, global_cols);
    return distribution;
}

void check_linearization(const int global_rows,
                         const int global_cols,
                         const StorageMajorOrder major_order)
{
    std::shared_ptr<BlockCyclicDistributionInformation> distribution
        = make_serial_distribution(global_rows, global_cols);
    const int local_extent
        = major_order == StorageMajorOrder::ColumnMajor ? global_rows : global_cols;
    const int ld = local_extent > 1 ? local_extent : 1;

    DistributedNDArray<double> array(distribution, major_order, ld);

    EXPECT_EQ(array.distribution().get(), distribution.get());
    EXPECT_EQ(array.major_order(), major_order);
    EXPECT_EQ(array.ld(), ld);
    EXPECT_EQ(array.nrloc(), global_rows);
    EXPECT_EQ(array.ncloc(), global_cols);
    EXPECT_EQ(array.local_block().data(), array.data());

    const std::size_t ld_size = static_cast<std::size_t>(array.ld());

    for (int local_row = 0; local_row < global_rows; ++local_row)
    {
        for (int local_col = 0; local_col < global_cols; ++local_col)
        {
            const double element_value = static_cast<double>(local_row * global_cols + local_col);
            array.at_local(local_row, local_col) = element_value;
        }
    }

    for (int local_row = 0; local_row < global_rows; ++local_row)
    {
        for (int local_col = 0; local_col < global_cols; ++local_col)
        {
            std::size_t expected_flat_index = 0;
            if (major_order == StorageMajorOrder::ColumnMajor)
            {
                expected_flat_index = static_cast<std::size_t>(local_col) * ld_size
                    + static_cast<std::size_t>(local_row);
            }
            else
            {
                expected_flat_index = static_cast<std::size_t>(local_row) * ld_size
                    + static_cast<std::size_t>(local_col);
            }

            const double* const expected_address = array.data() + expected_flat_index;
            EXPECT_EQ(&array.at_local(local_row, local_col), expected_address);

            const double element_value = static_cast<double>(local_row * global_cols + local_col);
            EXPECT_EQ(array.at_local(local_row, local_col), element_value);

            const int mapped_local_row = distribution->irloc(local_row);
            const int mapped_local_col = distribution->icloc(local_col);
            ASSERT_NE(mapped_local_row, -1);
            ASSERT_NE(mapped_local_col, -1);
            EXPECT_EQ(&array.at_global(local_row, local_col),
                      &array.at_local(mapped_local_row, mapped_local_col));
            EXPECT_EQ(array.at_global(local_row, local_col),
                      array.at_local(mapped_local_row, mapped_local_col));

            const bool expected_owner = distribution->in_this_processor(local_row, local_col);
            EXPECT_EQ(array.owns_global_index(local_row, local_col), expected_owner);
        }
    }
}
} // namespace

TEST(DistributedNDArrayTest, ColumnMajorLinearization)
{
    const int global_rows = 8;
    const int global_cols = 6;
    check_linearization(global_rows, global_cols, StorageMajorOrder::ColumnMajor);
}

TEST(DistributedNDArrayTest, RowMajorLinearization)
{
    const int global_rows = 2;
    const int global_cols = 3;
    check_linearization(global_rows, global_cols, StorageMajorOrder::RowMajor);
}

TEST(DistributedNDArrayTest, AtGlobalAssertsOutsideGlobalBounds)
{
    const int global_rows = 4;
    const int global_cols = 3;
    std::shared_ptr<BlockCyclicDistributionInformation> distribution
        = make_serial_distribution(global_rows, global_cols);
    const int ld = global_rows > 1 ? global_rows : 1;
    DistributedNDArray<double> array(distribution, StorageMajorOrder::ColumnMajor, ld);

    EXPECT_DEATH(array.at_global(global_rows, 0), "");
}

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
