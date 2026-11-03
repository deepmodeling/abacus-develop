#include "source_base/block_cyclic_distribution_information.h"
#include "source_base/parallel_common.h"

#include "gtest/gtest.h"

#include <algorithm>
#include <cstddef>
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
constexpr int DESCRIPTOR_ENTRY_COUNT = 9;
constexpr int DESCRIPTOR_MATRIX_TYPE_INDEX = 0;
constexpr int DESCRIPTOR_BLACS_CONTEXT_INDEX = 1;
constexpr int DESCRIPTOR_GLOBAL_ROWS_INDEX = 2;
constexpr int DESCRIPTOR_GLOBAL_COLS_INDEX = 3;
constexpr int DESCRIPTOR_BLOCK_ROWS_INDEX = 4;
constexpr int DESCRIPTOR_BLOCK_COLS_INDEX = 5;
constexpr int DESCRIPTOR_PROCESS_ROWS_INDEX = 6;
constexpr int DESCRIPTOR_PROCESS_COLS_INDEX = 7;
constexpr int DESCRIPTOR_LEADING_DIM_INDEX = 8;
constexpr int UNOWNED_BLACS_CONTEXT = -1;
} // namespace

TEST(BlockCyclicDistributionInformationTest, DefaultConstruction)
{
    BlockCyclicDistributionInformation distribution;

    EXPECT_EQ(distribution.nrloc(), 0);
    EXPECT_EQ(distribution.ncloc(), 0);
    EXPECT_EQ(distribution.size(), 0);
    EXPECT_EQ(distribution.nb(), 1);
    EXPECT_EQ(distribution.nrproc(), 0);
    EXPECT_EQ(distribution.ncproc(), 0);
    EXPECT_EQ(distribution.irproc(), UNOWNED_BLACS_CONTEXT);
    EXPECT_EQ(distribution.icproc(), UNOWNED_BLACS_CONTEXT);
    EXPECT_FALSE(distribution.is_serial());
#ifdef __MPI
    EXPECT_EQ(distribution.ictxt(), UNOWNED_BLACS_CONTEXT);
    const int* const descriptor = distribution.desc();
    for (int i = 0; i < DESCRIPTOR_ENTRY_COUNT; ++i)
    {
        EXPECT_EQ(descriptor[i], 0);
    }
#endif
}

TEST(BlockCyclicDistributionInformationTest, SerialInitializationAndIndexMaps)
{
    BlockCyclicDistributionInformation distribution;
    distribution.set_serial(7, 5);

    EXPECT_TRUE(distribution.is_serial());
    EXPECT_EQ(distribution.nrloc(), 7);
    EXPECT_EQ(distribution.ncloc(), 5);
    EXPECT_EQ(distribution.size(), 35);
    EXPECT_EQ(distribution.nb(), 1);
    EXPECT_EQ(distribution.nrproc(), 1);
    EXPECT_EQ(distribution.ncproc(), 1);
    EXPECT_EQ(distribution.irproc(), 0);
    EXPECT_EQ(distribution.icproc(), 0);

    for (int iat = 0; iat < 7; ++iat)
    {
        EXPECT_EQ(distribution.irloc(iat), iat);
        EXPECT_EQ(distribution.irglb(iat), iat);
    }
    for (int iat = 0; iat < 5; ++iat)
    {
        EXPECT_EQ(distribution.icloc(iat), iat);
        EXPECT_EQ(distribution.icglb(iat), iat);
    }
    EXPECT_TRUE(distribution.in_this_processor(6, 4));
    EXPECT_EQ(distribution.owner_processor(6, 4), 0);

    distribution.set_iproc(2, 3);
    EXPECT_EQ(distribution.irproc(), 2);
    EXPECT_EQ(distribution.icproc(), 3);
    EXPECT_TRUE(distribution.blacs_in_this_processor(2, 3));
    EXPECT_FALSE(distribution.blacs_in_this_processor(3, 2));
}

TEST(BlockCyclicDistributionInformationTest, SerialMovePreservesValues)
{
    BlockCyclicDistributionInformation source;
    source.set_serial(4, 6);

    BlockCyclicDistributionInformation moved(std::move(source));
    EXPECT_EQ(moved.nrloc(), 4);
    EXPECT_EQ(moved.ncloc(), 6);
    EXPECT_EQ(moved.size(), 24);
    EXPECT_TRUE(moved.is_serial());

    BlockCyclicDistributionInformation assigned;
    assigned = std::move(moved);
    EXPECT_EQ(assigned.nrloc(), 4);
    EXPECT_EQ(assigned.ncloc(), 6);
    EXPECT_EQ(assigned.size(), 24);
}

#ifdef __MPI
TEST(BlockCyclicDistributionInformationTest, MpiInitializationAndIndexMaps)
{
    BlockCyclicDistributionInformation distribution;
    const int global_rows = 30;
    const int global_cols = 35;
    const int block_size = 2;
    const int is_degenerate = distribution.init(global_rows, global_cols, block_size, MPI_COMM_WORLD, false);
    EXPECT_EQ(is_degenerate, 0);

    EXPECT_EQ(distribution.nrglb(), global_rows);
    EXPECT_EQ(distribution.ncglb(), global_cols);
    EXPECT_EQ(distribution.size(),
              static_cast<std::size_t>(distribution.nrloc())
                  * distribution.ncloc());
    const int process_count = Parallel_Common::comm_size(MPI_COMM_WORLD);
    EXPECT_EQ(distribution.nrproc() * distribution.ncproc(), process_count);

    const int* const descriptor = distribution.desc();
    ASSERT_NE(descriptor, nullptr);
    EXPECT_EQ(descriptor[DESCRIPTOR_MATRIX_TYPE_INDEX], 1);
    EXPECT_EQ(descriptor[DESCRIPTOR_BLACS_CONTEXT_INDEX], distribution.ictxt());
    EXPECT_EQ(descriptor[DESCRIPTOR_GLOBAL_ROWS_INDEX], global_rows);
    EXPECT_EQ(descriptor[DESCRIPTOR_GLOBAL_COLS_INDEX], global_cols);
    EXPECT_EQ(descriptor[DESCRIPTOR_BLOCK_ROWS_INDEX], block_size);
    EXPECT_EQ(descriptor[DESCRIPTOR_BLOCK_COLS_INDEX], block_size);
    EXPECT_EQ(descriptor[DESCRIPTOR_PROCESS_ROWS_INDEX], 0);
    EXPECT_EQ(descriptor[DESCRIPTOR_PROCESS_COLS_INDEX], 0);
    EXPECT_EQ(descriptor[DESCRIPTOR_LEADING_DIM_INDEX],
              std::max(distribution.nrloc(), 1));
    for (int iat = 0; iat < global_rows; ++iat)
    {
        const int local_row = distribution.irloc(iat);
        if (local_row != -1)
        {
            EXPECT_EQ(distribution.irglb(local_row), iat);
        }
    }
    for (int iat = 0; iat < global_cols; ++iat)
    {
        const int local_col = distribution.icloc(iat);
        if (local_col != -1)
        {
            EXPECT_EQ(distribution.icglb(local_col), iat);
        }
    }

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    EXPECT_EQ(distribution.in_this_processor(1, 1),
              distribution.owner_processor(1, 1) == rank);
    EXPECT_TRUE(distribution.blacs_in_this_processor(distribution.irproc(),
                                              distribution.icproc()));
}

TEST(BlockCyclicDistributionInformationTest, MpiMoveTransfersBlacsOwnershipOnce)
{
    BlockCyclicDistributionInformation source;
    source.init(12, 10, 2, MPI_COMM_WORLD, false);
    const int source_context = source.ictxt();
    ASSERT_GE(source_context, 0);

    BlockCyclicDistributionInformation moved(std::move(source));
    EXPECT_EQ(moved.ictxt(), source_context);
    EXPECT_EQ(source.ictxt(), UNOWNED_BLACS_CONTEXT);
}

TEST(BlockCyclicDistributionInformationTest, BindDescriptorBorrowsBlacsContext)
{
    BlockCyclicDistributionInformation owner;
    owner.init(12, 10, 2, MPI_COMM_WORLD, false);
    const int owner_context = owner.ictxt();
    ASSERT_GE(owner_context, 0);

    BlockCyclicDistributionInformation borrower;
    borrower.bind_descriptor(8, 8, 2, owner_context);
    EXPECT_EQ(borrower.ictxt(), owner_context);

    BlockCyclicDistributionInformation moved_borrower(std::move(borrower));
    EXPECT_EQ(moved_borrower.ictxt(), owner_context);
    EXPECT_EQ(borrower.ictxt(), UNOWNED_BLACS_CONTEXT);
    EXPECT_GE(owner.ictxt(), 0);
}

TEST(BlockCyclicDistributionInformationTest, SharedBorrowedContextSupportsMultipleLayouts)
{
    BlockCyclicDistributionInformation owner;
    owner.init(30, 35, 2, MPI_COMM_WORLD, false);
    const int owner_context = owner.ictxt();
    ASSERT_GE(owner_context, 0);

    BlockCyclicDistributionInformation first_borrower;
    const int first_degenerate = first_borrower.bind_descriptor(30, 35, 2, owner_context);
    EXPECT_EQ(first_degenerate, 0);

    BlockCyclicDistributionInformation second_borrower;
    const int second_degenerate = second_borrower.bind_descriptor(49, 94, 3, owner_context);
    EXPECT_EQ(second_degenerate, 0);

    EXPECT_EQ(first_borrower.ictxt(), owner_context);
    EXPECT_EQ(second_borrower.ictxt(), owner_context);
    EXPECT_EQ(first_borrower.desc()[DESCRIPTOR_BLACS_CONTEXT_INDEX], owner_context);
    EXPECT_EQ(second_borrower.desc()[DESCRIPTOR_BLACS_CONTEXT_INDEX], owner_context);

    BlockCyclicDistributionInformation independent_owner;
    independent_owner.init(57, 57, 2, MPI_COMM_WORLD, false);
    EXPECT_NE(independent_owner.ictxt(), owner_context);
}

TEST(BlockCyclicDistributionInformationTest, RebindingOwnedContextKeepsGridAlive)
{
    BlockCyclicDistributionInformation distribution;
    distribution.init(30, 35, 2, MPI_COMM_WORLD, false);
    const int original_context = distribution.ictxt();
    ASSERT_GE(original_context, 0);

    const int is_degenerate = distribution.bind_descriptor(30, 35, 1, original_context);
    EXPECT_EQ(is_degenerate, 0);

    const int* const descriptor = distribution.desc();
    EXPECT_EQ(distribution.ictxt(), original_context);
    EXPECT_EQ(descriptor[DESCRIPTOR_MATRIX_TYPE_INDEX], 1);
    EXPECT_EQ(descriptor[DESCRIPTOR_BLACS_CONTEXT_INDEX], original_context);
    EXPECT_EQ(descriptor[DESCRIPTOR_GLOBAL_ROWS_INDEX], 30);
    EXPECT_EQ(descriptor[DESCRIPTOR_GLOBAL_COLS_INDEX], 35);
    EXPECT_EQ(descriptor[DESCRIPTOR_BLOCK_ROWS_INDEX], 1);
    EXPECT_EQ(descriptor[DESCRIPTOR_BLOCK_COLS_INDEX], 1);
    EXPECT_EQ(distribution.nb(), 1);
    EXPECT_EQ(distribution.size(),
              static_cast<std::size_t>(distribution.nrloc())
                  * distribution.ncloc());
    EXPECT_EQ(distribution.nrproc() * distribution.ncproc(),
              Parallel_Common::comm_size(MPI_COMM_WORLD));
}

TEST(BlockCyclicDistributionInformationTest, SerialInitializationClearsBorrowedContext)
{
    BlockCyclicDistributionInformation owner;
    owner.init(30, 35, 1, MPI_COMM_WORLD, false);
    const int owner_context = owner.ictxt();
    ASSERT_GE(owner_context, 0);

    BlockCyclicDistributionInformation borrower;
    const int is_degenerate = borrower.bind_descriptor(30, 35, 1, owner_context);
    EXPECT_EQ(is_degenerate, 0);
    EXPECT_EQ(borrower.ictxt(), owner_context);

    borrower.set_serial(3, 4);
    EXPECT_EQ(borrower.ictxt(), UNOWNED_BLACS_CONTEXT);
    EXPECT_EQ(borrower.comm(), MPI_COMM_NULL);
    EXPECT_GE(owner.ictxt(), 0);
    EXPECT_NE(owner.comm(), MPI_COMM_NULL);
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
