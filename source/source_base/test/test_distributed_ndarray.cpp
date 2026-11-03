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
#ifdef __MPI
constexpr int DESCRIPTOR_MATRIX_TYPE_INDEX = 0;
constexpr int DESCRIPTOR_BLACS_CONTEXT_INDEX = 1;
constexpr int DESCRIPTOR_GLOBAL_ROWS_INDEX = 2;
constexpr int DESCRIPTOR_GLOBAL_COLS_INDEX = 3;
constexpr int DESCRIPTOR_BLOCK_ROWS_INDEX = 4;
constexpr int DESCRIPTOR_BLOCK_COLS_INDEX = 5;
constexpr int DESCRIPTOR_PROCESS_ROWS_INDEX = 6;
constexpr int DESCRIPTOR_PROCESS_COLS_INDEX = 7;
constexpr int DESCRIPTOR_LEADING_DIM_INDEX = 8;
#endif

std::shared_ptr<BlockCyclicDistributionInformation> make_serial_distribution(
    const int global_rows,
    const int global_cols)
{
    std::shared_ptr<BlockCyclicDistributionInformation> distribution
        = std::make_shared<BlockCyclicDistributionInformation>();
    distribution->set_serial(global_rows, global_cols);

#ifdef __MPI
    int* const descriptor = distribution->desc();
    descriptor[DESCRIPTOR_MATRIX_TYPE_INDEX] = 1;
    descriptor[DESCRIPTOR_BLACS_CONTEXT_INDEX] = -1;
    descriptor[DESCRIPTOR_GLOBAL_ROWS_INDEX] = global_rows;
    descriptor[DESCRIPTOR_GLOBAL_COLS_INDEX] = global_cols;
    descriptor[DESCRIPTOR_BLOCK_ROWS_INDEX] = 1;
    descriptor[DESCRIPTOR_BLOCK_COLS_INDEX] = 1;
    descriptor[DESCRIPTOR_PROCESS_ROWS_INDEX] = 0;
    descriptor[DESCRIPTOR_PROCESS_COLS_INDEX] = 0;
    descriptor[DESCRIPTOR_LEADING_DIM_INDEX] = global_rows > 1 ? global_rows : 1;
#endif

    return distribution;
}

void check_linearization(const int global_rows,
                         const int global_cols,
                         const StorageMajorOrder major_order)
{
    std::shared_ptr<BlockCyclicDistributionInformation> distribution
        = make_serial_distribution(global_rows, global_cols);
    const int expected_leading_dimension = global_rows > 1 ? global_rows : 1;

    DistributedNDArray<double> array(distribution, major_order, expected_leading_dimension);

    EXPECT_EQ(array.distribution().get(), distribution.get());
    EXPECT_EQ(array.major_order(), major_order);
    EXPECT_EQ(array.leading_dimension(), expected_leading_dimension);
    EXPECT_EQ(array.nrloc(), global_rows);
    EXPECT_EQ(array.ncloc(), global_cols);
    EXPECT_EQ(array.local_block().data(), array.data());

    std::size_t leading_dimension = static_cast<std::size_t>(array.leading_dimension());
#ifdef __MPI
    const int* const descriptor = distribution->desc();
    EXPECT_EQ(array.leading_dimension(), descriptor[DESCRIPTOR_LEADING_DIM_INDEX]);
    leading_dimension = static_cast<std::size_t>(descriptor[DESCRIPTOR_LEADING_DIM_INDEX]);
#endif

    for (int local_row = 0; local_row < global_rows; ++local_row)
    {
        for (int local_col = 0; local_col < global_cols; ++local_col)
        {
            std::size_t expected_flat_index = 0;
            if (major_order == StorageMajorOrder::ColumnMajor)
            {
                expected_flat_index = static_cast<std::size_t>(local_col) * leading_dimension
                    + static_cast<std::size_t>(local_row);
            }
            else
            {
                expected_flat_index = static_cast<std::size_t>(local_row) * leading_dimension
                    + static_cast<std::size_t>(local_col);
            }

            const double* const expected_address = array.data() + expected_flat_index;
            EXPECT_EQ(&array.at_local(local_row, local_col), expected_address);

            const double element_value = static_cast<double>(local_row * global_cols + local_col);
            array.at_local(local_row, local_col) = element_value;
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
    const int global_rows = 5;
    const int global_cols = 3;
    check_linearization(global_rows, global_cols, StorageMajorOrder::RowMajor);
}

TEST(DistributedNDArrayTest, AtGlobalAssertsOutsideGlobalBounds)
{
    const int global_rows = 4;
    const int global_cols = 3;
    std::shared_ptr<BlockCyclicDistributionInformation> distribution
        = make_serial_distribution(global_rows, global_cols);
    const int leading_dimension = global_rows > 1 ? global_rows : 1;
    DistributedNDArray<double> array(distribution, StorageMajorOrder::ColumnMajor, leading_dimension);

    EXPECT_DEATH(array.at_global(global_rows, 0), "");
}

#ifdef __MPI
TEST(DistributedNDArrayTest, RejectsMismatchedDescriptorLeadingDimension)
{
    const int global_rows = 4;
    const int global_cols = 3;
    std::shared_ptr<BlockCyclicDistributionInformation> distribution
        = make_serial_distribution(global_rows, global_cols);
    int* const descriptor = distribution->desc();
    descriptor[DESCRIPTOR_LEADING_DIM_INDEX] += 1;
    const int expected_leading_dimension = global_rows > 1 ? global_rows : 1;

    EXPECT_DEATH(
        {
            DistributedNDArray<double> invalid_array(
                distribution,
                StorageMajorOrder::ColumnMajor,
                expected_leading_dimension);
            static_cast<void>(invalid_array);
        },
        "");
}
#endif

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
