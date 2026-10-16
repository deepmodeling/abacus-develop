#include <gtest/gtest.h>

#include "source_cell/module_neighlist/page_allocator.h"

#include <stdexcept>

TEST(PageAllocatorTest, RejectsNonPositivePageSize)
{
    PageAllocator allocator;

    EXPECT_THROW(allocator.initialize(0), std::invalid_argument);
    EXPECT_THROW(allocator.initialize(-1), std::invalid_argument);
}

TEST(PageAllocatorTest, AllocatesContiguousBlocksWithinPage)
{
    PageAllocator allocator;
    allocator.initialize(8);

    int* first = allocator.allocate(2);
    int* second = allocator.allocate(3);

    ASSERT_NE(first, nullptr);
    ASSERT_NE(second, nullptr);
    EXPECT_EQ(second, first + 2);
    EXPECT_EQ(allocator.get_pgsize(), 8);
}

TEST(PageAllocatorTest, StartsNewPageAndResetReusesFirstPage)
{
    PageAllocator allocator;
    allocator.initialize(2);

    int* first = allocator.allocate(2);
    int* second_page = allocator.allocate(1);

    ASSERT_NE(first, nullptr);
    ASSERT_NE(second_page, nullptr);
    EXPECT_NE(second_page, first);

    allocator.reset();
    EXPECT_EQ(allocator.allocate(1), first);
}

TEST(PageAllocatorTest, NonPositiveAllocationReturnsNull)
{
    PageAllocator allocator;
    allocator.initialize(4);

    EXPECT_EQ(allocator.allocate(0), nullptr);
    EXPECT_EQ(allocator.allocate(-1), nullptr);
}
