#include "source_base/parallel_blkcyc.h"
#include "source_base/parallel_common.h"

#include "gtest/gtest.h"

#include <algorithm>
#include <cstddef>
#include <type_traits>
#include <utility>
#include <vector>

static_assert(std::is_default_constructible<BlockCyclicDistribution>::value,
              "BlockCyclicDistribution must be default constructible");
static_assert(!std::is_copy_constructible<BlockCyclicDistribution>::value,
              "BlockCyclicDistribution must not be copy constructible");
static_assert(!std::is_copy_assignable<BlockCyclicDistribution>::value,
              "BlockCyclicDistribution must not be copy assignable");
static_assert(std::is_move_constructible<BlockCyclicDistribution>::value,
              "BlockCyclicDistribution must be move constructible");
static_assert(std::is_move_assignable<BlockCyclicDistribution>::value,
              "BlockCyclicDistribution must be move assignable");
static_assert(std::is_nothrow_move_constructible<BlockCyclicDistribution>::value,
              "BlockCyclicDistribution move construction must be noexcept");
static_assert(std::is_nothrow_move_assignable<BlockCyclicDistribution>::value,
              "BlockCyclicDistribution move assignment must be noexcept");

namespace
{
constexpr int DESC_SIZE = 9;
constexpr int DESC_DTYPE = 0;
constexpr int DESC_ICTXT = 1;
constexpr int DESC_M = 2;
constexpr int DESC_N = 3;
constexpr int DESC_MB = 4;
constexpr int DESC_NB = 5;
constexpr int DESC_RSRC = 6;
constexpr int DESC_CSRC = 7;
constexpr int DESC_LLD = 8;
constexpr int UNOWNED_BLACS_CONTEXT = -1;
} // namespace

TEST(BlockCyclicDistributionTest, DefaultConstruction)
{
    BlockCyclicDistribution distribution;

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
    for (int i = 0; i < DESC_SIZE; ++i)
    {
        EXPECT_EQ(descriptor[i], 0);
    }
#endif
}

TEST(BlockCyclicDistributionTest, SerialInitializationAndIndexMaps)
{
    BlockCyclicDistribution distribution;
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

TEST(BlockCyclicDistributionTest, SerialMovePreservesValues)
{
    BlockCyclicDistribution source;
    source.set_serial(4, 6);

    BlockCyclicDistribution moved(std::move(source));
    EXPECT_EQ(moved.nrloc(), 4);
    EXPECT_EQ(moved.ncloc(), 6);
    EXPECT_EQ(moved.size(), 24);
    EXPECT_TRUE(moved.is_serial());

    BlockCyclicDistribution assigned;
    assigned = std::move(moved);
    EXPECT_EQ(assigned.nrloc(), 4);
    EXPECT_EQ(assigned.ncloc(), 6);
    EXPECT_EQ(assigned.size(), 24);
}

#ifdef __MPI
TEST(BlockCyclicDistributionTest, MpiInitializationAndIndexMaps)
{
    BlockCyclicDistribution distribution;
    const int nrglb = 30;
    const int ncglb = 35;
    const int block_size = 2;
    const int is_degenerate = distribution.init(nrglb, ncglb, block_size, MPI_COMM_WORLD, false);
    EXPECT_EQ(is_degenerate, 0);

    EXPECT_EQ(distribution.nrglb(), nrglb);
    EXPECT_EQ(distribution.ncglb(), ncglb);
    EXPECT_EQ(distribution.size(),
              static_cast<std::size_t>(distribution.nrloc())
                  * distribution.ncloc());
    const int process_count = Parallel_Common::comm_size(MPI_COMM_WORLD);
    EXPECT_EQ(distribution.nrproc() * distribution.ncproc(), process_count);

    const int* const descriptor = distribution.desc();
    ASSERT_NE(descriptor, nullptr);
    EXPECT_EQ(descriptor[DESC_DTYPE], 1);
    EXPECT_EQ(descriptor[DESC_ICTXT], distribution.ictxt());
    EXPECT_EQ(descriptor[DESC_M], nrglb);
    EXPECT_EQ(descriptor[DESC_N], ncglb);
    EXPECT_EQ(descriptor[DESC_MB], block_size);
    EXPECT_EQ(descriptor[DESC_NB], block_size);
    EXPECT_EQ(descriptor[DESC_RSRC], 0);
    EXPECT_EQ(descriptor[DESC_CSRC], 0);
    EXPECT_EQ(descriptor[DESC_LLD],
              std::max(distribution.nrloc(), 1));
    for (int i = 0; i < nrglb; ++i)
    {
        const int il = distribution.irloc(i);
        if (il != -1)
        {
            EXPECT_EQ(distribution.irglb(il), i);
        }
    }
    for (int j = 0; j < ncglb; ++j)
    {
        const int jl = distribution.icloc(j);
        if (jl != -1)
        {
            EXPECT_EQ(distribution.icglb(jl), j);
        }
    }

    int irank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &irank);
    EXPECT_EQ(distribution.in_this_processor(1, 1),
              distribution.owner_processor(1, 1) == irank);
    EXPECT_TRUE(distribution.blacs_in_this_processor(distribution.irproc(),
                                              distribution.icproc()));
}

TEST(BlockCyclicDistributionTest, MpiMoveTransfersBlacsOwnershipOnce)
{
    BlockCyclicDistribution source;
    source.init(12, 10, 2, MPI_COMM_WORLD, false);
    const int source_context = source.ictxt();
    ASSERT_GE(source_context, 0);

    BlockCyclicDistribution moved(std::move(source));
    EXPECT_EQ(moved.ictxt(), source_context);
    EXPECT_EQ(source.ictxt(), UNOWNED_BLACS_CONTEXT);
}

TEST(BlockCyclicDistributionTest, BindDescriptorBorrowsBlacsContext)
{
    BlockCyclicDistribution owner;
    owner.init(12, 10, 2, MPI_COMM_WORLD, false);
    const int owner_context = owner.ictxt();
    ASSERT_GE(owner_context, 0);

    BlockCyclicDistribution borrower;
    borrower.bind_descriptor(8, 8, 2, owner_context);
    EXPECT_EQ(borrower.ictxt(), owner_context);

    BlockCyclicDistribution moved_borrower(std::move(borrower));
    EXPECT_EQ(moved_borrower.ictxt(), owner_context);
    EXPECT_EQ(borrower.ictxt(), UNOWNED_BLACS_CONTEXT);
    EXPECT_GE(owner.ictxt(), 0);
}

TEST(BlockCyclicDistributionTest, SharedBorrowedContextSupportsMultipleLayouts)
{
    BlockCyclicDistribution owner;
    owner.init(30, 35, 2, MPI_COMM_WORLD, false);
    const int owner_context = owner.ictxt();
    ASSERT_GE(owner_context, 0);

    BlockCyclicDistribution first_borrower;
    const int first_degenerate = first_borrower.bind_descriptor(30, 35, 2, owner_context);
    EXPECT_EQ(first_degenerate, 0);

    BlockCyclicDistribution second_borrower;
    const int second_degenerate = second_borrower.bind_descriptor(49, 94, 3, owner_context);
    EXPECT_EQ(second_degenerate, 0);

    EXPECT_EQ(first_borrower.ictxt(), owner_context);
    EXPECT_EQ(second_borrower.ictxt(), owner_context);
    EXPECT_EQ(first_borrower.desc()[DESC_ICTXT], owner_context);
    EXPECT_EQ(second_borrower.desc()[DESC_ICTXT], owner_context);

    BlockCyclicDistribution independent_owner;
    independent_owner.init(57, 57, 2, MPI_COMM_WORLD, false);
    EXPECT_NE(independent_owner.ictxt(), owner_context);
}

TEST(BlockCyclicDistributionTest, RebindingOwnedContextKeepsGridAlive)
{
    BlockCyclicDistribution distribution;
    distribution.init(30, 35, 2, MPI_COMM_WORLD, false);
    const int original_context = distribution.ictxt();
    ASSERT_GE(original_context, 0);

    const int is_degenerate = distribution.bind_descriptor(30, 35, 1, original_context);
    EXPECT_EQ(is_degenerate, 0);

    const int* const descriptor = distribution.desc();
    EXPECT_EQ(distribution.ictxt(), original_context);
    EXPECT_EQ(descriptor[DESC_DTYPE], 1);
    EXPECT_EQ(descriptor[DESC_ICTXT], original_context);
    EXPECT_EQ(descriptor[DESC_M], 30);
    EXPECT_EQ(descriptor[DESC_N], 35);
    EXPECT_EQ(descriptor[DESC_MB], 1);
    EXPECT_EQ(descriptor[DESC_NB], 1);
    EXPECT_EQ(distribution.nb(), 1);
    EXPECT_EQ(distribution.size(),
              static_cast<std::size_t>(distribution.nrloc())
                  * distribution.ncloc());
    EXPECT_EQ(distribution.nrproc() * distribution.ncproc(),
              Parallel_Common::comm_size(MPI_COMM_WORLD));
}

TEST(BlockCyclicDistributionTest, SerialInitializationClearsBorrowedContext)
{
    BlockCyclicDistribution owner;
    owner.init(30, 35, 1, MPI_COMM_WORLD, false);
    const int owner_context = owner.ictxt();
    ASSERT_GE(owner_context, 0);

    BlockCyclicDistribution borrower;
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
