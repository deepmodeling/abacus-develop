#include "source_base/block_cyclic_distribution_information.h"
#include "source_base/parallel_common.h"

#include "gtest/gtest.h"

#include <algorithm>
#include <cstdint>
#include <type_traits>
#include <utility>
#include <vector>

static_assert(std::is_default_constructible<BlockCyclicDistributionInformation>::value,
              "BlockCyclicDistributionInformation must be default constructible");
static_assert(!std::is_copy_constructible<BlockCyclicDistributionInformation>::value,
              "BlockCyclicDistributionInformation must not be copy constructible");
static_assert(!std::is_copy_assignable<BlockCyclicDistributionInformation>::value,
              "BlockCyclicDistributionInformation must not be copy assignable");
static_assert(std::is_move_constructible<BlockCyclicDistributionInformation>::value,
              "BlockCyclicDistributionInformation must be move constructible");
static_assert(std::is_move_assignable<BlockCyclicDistributionInformation>::value,
              "BlockCyclicDistributionInformation must be move assignable");
static_assert(std::is_nothrow_move_constructible<BlockCyclicDistributionInformation>::value,
              "BlockCyclicDistributionInformation move construction must be noexcept");
static_assert(std::is_nothrow_move_assignable<BlockCyclicDistributionInformation>::value,
              "BlockCyclicDistributionInformation move assignment must be noexcept");

namespace
{
constexpr int kDescriptorEntryCount = 9;
constexpr int kDescriptorMatrixTypeIndex = 0;
constexpr int kDescriptorBlacsContextIndex = 1;
constexpr int kDescriptorGlobalRowsIndex = 2;
constexpr int kDescriptorGlobalColsIndex = 3;
constexpr int kDescriptorBlockRowsIndex = 4;
constexpr int kDescriptorBlockColsIndex = 5;
constexpr int kDescriptorProcessRowsIndex = 6;
constexpr int kDescriptorProcessColsIndex = 7;
constexpr int kDescriptorLeadingDimIndex = 8;
constexpr int kUnownedBlacsContext = -1;
} // namespace

TEST(BlockCyclicDistributionInformationTest, DefaultConstruction)
{
    BlockCyclicDistributionInformation distribution;

    EXPECT_EQ(distribution.get_number_of_local_rows(), 0);
    EXPECT_EQ(distribution.get_number_of_local_cols(), 0);
    EXPECT_EQ(distribution.get_number_of_local_elements(), 0);
    EXPECT_EQ(distribution.get_block_size(), 1);
    EXPECT_EQ(distribution.get_grid_rows(), 0);
    EXPECT_EQ(distribution.get_grid_cols(), 0);
    EXPECT_EQ(distribution.get_coord_row(), kUnownedBlacsContext);
    EXPECT_EQ(distribution.get_coord_col(), kUnownedBlacsContext);
    EXPECT_FALSE(distribution.is_serial());
#ifdef __MPI
    EXPECT_EQ(distribution.get_blacs_context(), kUnownedBlacsContext);
    const int* const descriptor = distribution.get_descriptor();
    for (int i = 0; i < kDescriptorEntryCount; ++i)
    {
        EXPECT_EQ(descriptor[i], 0);
    }
#endif
}

TEST(BlockCyclicDistributionInformationTest, SerialInitializationAndIndexMaps)
{
    BlockCyclicDistributionInformation distribution;
    distribution.initialize_serial(7, 5);

    EXPECT_TRUE(distribution.is_serial());
    EXPECT_EQ(distribution.get_number_of_local_rows(), 7);
    EXPECT_EQ(distribution.get_number_of_local_cols(), 5);
    EXPECT_EQ(distribution.get_number_of_local_elements(), 35);
    EXPECT_EQ(distribution.get_block_size(), 1);
    EXPECT_EQ(distribution.get_grid_rows(), 1);
    EXPECT_EQ(distribution.get_grid_cols(), 1);
    EXPECT_EQ(distribution.get_coord_row(), 0);
    EXPECT_EQ(distribution.get_coord_col(), 0);

    for (int iat = 0; iat < 7; ++iat)
    {
        EXPECT_EQ(distribution.global_to_local_row(iat), iat);
        EXPECT_EQ(distribution.local_to_global_row(iat), iat);
    }
    for (int iat = 0; iat < 5; ++iat)
    {
        EXPECT_EQ(distribution.global_to_local_col(iat), iat);
        EXPECT_EQ(distribution.local_to_global_col(iat), iat);
    }
    EXPECT_TRUE(distribution.owns_index(6, 4));
    EXPECT_EQ(distribution.owner_of_index(6, 4), 0);

    distribution.set_coord(2, 3);
    EXPECT_EQ(distribution.get_coord_row(), 2);
    EXPECT_EQ(distribution.get_coord_col(), 3);
    EXPECT_TRUE(distribution.owns_blacs_index(2, 3));
    EXPECT_FALSE(distribution.owns_blacs_index(3, 2));
}

TEST(BlockCyclicDistributionInformationTest, SerialMovePreservesValues)
{
    BlockCyclicDistributionInformation source;
    source.initialize_serial(4, 6);

    BlockCyclicDistributionInformation moved(std::move(source));
    EXPECT_EQ(moved.get_number_of_local_rows(), 4);
    EXPECT_EQ(moved.get_number_of_local_cols(), 6);
    EXPECT_EQ(moved.get_number_of_local_elements(), 24);
    EXPECT_TRUE(moved.is_serial());

    BlockCyclicDistributionInformation assigned;
    assigned = std::move(moved);
    EXPECT_EQ(assigned.get_number_of_local_rows(), 4);
    EXPECT_EQ(assigned.get_number_of_local_cols(), 6);
    EXPECT_EQ(assigned.get_number_of_local_elements(), 24);
}

#ifdef __MPI
TEST(BlockCyclicDistributionInformationTest, MpiInitializationAndIndexMaps)
{
    BlockCyclicDistributionInformation distribution;
    const int global_rows = 30;
    const int global_cols = 35;
    const int block_size = 2;
    const int is_degenerate = distribution.initialize(global_rows, global_cols, block_size, MPI_COMM_WORLD);
    EXPECT_EQ(is_degenerate, 0);

    EXPECT_EQ(distribution.get_number_of_global_rows(), global_rows);
    EXPECT_EQ(distribution.get_number_of_global_cols(), global_cols);
    EXPECT_EQ(distribution.get_number_of_local_elements(),
              static_cast<int64_t>(distribution.get_number_of_local_rows())
                  * distribution.get_number_of_local_cols());
    const int process_count = Parallel_Common::comm_size(MPI_COMM_WORLD);
    EXPECT_EQ(distribution.get_grid_rows() * distribution.get_grid_cols(), process_count);

    const int* const descriptor = distribution.get_descriptor();
    ASSERT_NE(descriptor, nullptr);
    EXPECT_EQ(descriptor, distribution.get_mutable_descriptor());
    EXPECT_EQ(descriptor[kDescriptorMatrixTypeIndex], 1);
    EXPECT_EQ(descriptor[kDescriptorBlacsContextIndex], distribution.get_blacs_context());
    EXPECT_EQ(descriptor[kDescriptorGlobalRowsIndex], global_rows);
    EXPECT_EQ(descriptor[kDescriptorGlobalColsIndex], global_cols);
    EXPECT_EQ(descriptor[kDescriptorBlockRowsIndex], block_size);
    EXPECT_EQ(descriptor[kDescriptorBlockColsIndex], block_size);
    EXPECT_EQ(descriptor[kDescriptorProcessRowsIndex], 0);
    EXPECT_EQ(descriptor[kDescriptorProcessColsIndex], 0);
    EXPECT_EQ(descriptor[kDescriptorLeadingDimIndex],
              std::max(distribution.get_number_of_local_rows(), 1));
    for (int iat = 0; iat < global_rows; ++iat)
    {
        const int local_row = distribution.global_to_local_row(iat);
        if (local_row != -1)
        {
            EXPECT_EQ(distribution.local_to_global_row(local_row), iat);
        }
    }
    for (int iat = 0; iat < global_cols; ++iat)
    {
        const int local_col = distribution.global_to_local_col(iat);
        if (local_col != -1)
        {
            EXPECT_EQ(distribution.local_to_global_col(local_col), iat);
        }
    }

    EXPECT_EQ(distribution.owns_index(1, 1), distribution.in_this_processor(1, 1));
    EXPECT_EQ(distribution.owner_of_index(1, 1), distribution.owner_processor(1, 1));
    EXPECT_TRUE(distribution.owns_blacs_index(distribution.get_coord_row(),
                                              distribution.get_coord_col()));
}

TEST(BlockCyclicDistributionInformationTest, MpiMoveTransfersBlacsOwnershipOnce)
{
    BlockCyclicDistributionInformation source;
    source.initialize(12, 10, 2, MPI_COMM_WORLD);
    const int source_context = source.get_blacs_context();
    ASSERT_GE(source_context, 0);

    BlockCyclicDistributionInformation moved(std::move(source));
    EXPECT_EQ(moved.get_blacs_context(), source_context);
    EXPECT_EQ(source.get_blacs_context(), kUnownedBlacsContext);
}

TEST(BlockCyclicDistributionInformationTest, BindDescriptorBorrowsBlacsContext)
{
    BlockCyclicDistributionInformation owner;
    owner.initialize(12, 10, 2, MPI_COMM_WORLD);
    const int owner_context = owner.get_blacs_context();
    ASSERT_GE(owner_context, 0);

    BlockCyclicDistributionInformation borrower;
    borrower.bind_descriptor(8, 8, 2, owner_context);
    EXPECT_EQ(borrower.get_blacs_context(), owner_context);

    BlockCyclicDistributionInformation moved_borrower(std::move(borrower));
    EXPECT_EQ(moved_borrower.get_blacs_context(), owner_context);
    EXPECT_EQ(borrower.get_blacs_context(), kUnownedBlacsContext);
    EXPECT_GE(owner.get_blacs_context(), 0);
}

TEST(BlockCyclicDistributionInformationTest, SharedBorrowedContextSupportsMultipleLayouts)
{
    BlockCyclicDistributionInformation owner;
    owner.initialize(30, 35, 2, MPI_COMM_WORLD);
    const int owner_context = owner.get_blacs_context();
    ASSERT_GE(owner_context, 0);

    BlockCyclicDistributionInformation first_borrower;
    const int first_degenerate = first_borrower.bind_descriptor(30, 35, 2, owner_context);
    EXPECT_EQ(first_degenerate, 0);

    BlockCyclicDistributionInformation second_borrower;
    const int second_degenerate = second_borrower.bind_descriptor(49, 94, 3, owner_context);
    EXPECT_EQ(second_degenerate, 0);

    EXPECT_EQ(first_borrower.get_blacs_context(), owner_context);
    EXPECT_EQ(second_borrower.get_blacs_context(), owner_context);
    EXPECT_EQ(first_borrower.get_descriptor()[kDescriptorBlacsContextIndex], owner_context);
    EXPECT_EQ(second_borrower.get_descriptor()[kDescriptorBlacsContextIndex], owner_context);

    BlockCyclicDistributionInformation independent_owner;
    independent_owner.initialize(57, 57, 2, MPI_COMM_WORLD);
    EXPECT_NE(independent_owner.get_blacs_context(), owner_context);
}

TEST(BlockCyclicDistributionInformationTest, RebindingOwnedContextKeepsGridAlive)
{
    BlockCyclicDistributionInformation distribution;
    distribution.initialize(30, 35, 2, MPI_COMM_WORLD);
    const int original_context = distribution.get_blacs_context();
    ASSERT_GE(original_context, 0);

    const int is_degenerate = distribution.bind_descriptor(30, 35, 1, original_context);
    EXPECT_EQ(is_degenerate, 0);

    const int* const descriptor = distribution.get_descriptor();
    EXPECT_EQ(distribution.get_blacs_context(), original_context);
    EXPECT_EQ(descriptor[kDescriptorMatrixTypeIndex], 1);
    EXPECT_EQ(descriptor[kDescriptorBlacsContextIndex], original_context);
    EXPECT_EQ(descriptor[kDescriptorGlobalRowsIndex], 30);
    EXPECT_EQ(descriptor[kDescriptorGlobalColsIndex], 35);
    EXPECT_EQ(descriptor[kDescriptorBlockRowsIndex], 1);
    EXPECT_EQ(descriptor[kDescriptorBlockColsIndex], 1);
    EXPECT_EQ(distribution.get_block_size(), 1);
    EXPECT_EQ(distribution.get_number_of_local_elements(),
              static_cast<int64_t>(distribution.get_number_of_local_rows())
                  * distribution.get_number_of_local_cols());
    EXPECT_EQ(distribution.get_grid_rows() * distribution.get_grid_cols(),
              Parallel_Common::comm_size(MPI_COMM_WORLD));
}

TEST(BlockCyclicDistributionInformationTest, SerialInitializationClearsBorrowedContext)
{
    BlockCyclicDistributionInformation owner;
    owner.initialize(30, 35, 1, MPI_COMM_WORLD);
    const int owner_context = owner.get_blacs_context();
    ASSERT_GE(owner_context, 0);

    BlockCyclicDistributionInformation borrower;
    const int is_degenerate = borrower.bind_descriptor(30, 35, 1, owner_context);
    EXPECT_EQ(is_degenerate, 0);
    EXPECT_EQ(borrower.get_blacs_context(), owner_context);

    borrower.initialize_serial(3, 4);
    EXPECT_EQ(borrower.get_blacs_context(), kUnownedBlacsContext);
    EXPECT_EQ(borrower.get_mpi_communicator(), MPI_COMM_NULL);
    EXPECT_GE(owner.get_blacs_context(), 0);
    EXPECT_NE(owner.get_mpi_communicator(), MPI_COMM_NULL);
}
#endif

int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
#endif
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    MPI_Finalize();
#endif
    return result;
}
