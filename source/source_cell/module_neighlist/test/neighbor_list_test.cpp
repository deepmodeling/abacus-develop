#include <gtest/gtest.h>
#include "../neighbor_list.h"

// Full coverage tests for PageAllocator and NeighborList

TEST(PageAllocator_ConstructorsAndNewPage, DefaultAndCustom)
{
    PageAllocator pa_def; // default page size
    ASSERT_EQ(pa_def.pages.size(), 1u);
    EXPECT_EQ(pa_def.pages[0].capacity, pa_def.pgsize);

    PageAllocator pa(4);
    EXPECT_EQ(pa.pgsize, 4);
    size_t before = pa.pages.size();
    pa.new_page();
    EXPECT_EQ(pa.pages.size(), before + 1);
}

TEST(PageAllocator_AllocateEdgeCases, ZeroNegativeAndTooLarge)
{
    PageAllocator pa(8);
    EXPECT_EQ(pa.allocate(0), nullptr);
    EXPECT_EQ(pa.allocate(-5), nullptr);
    // larger than single page -> rejected
    EXPECT_EQ(pa.allocate(9), nullptr);
}

TEST(PageAllocator_AllocationBehavior, ExactPageAndOverflow)
{
    PageAllocator pa(4);
    int* p1 = pa.allocate(4);
    ASSERT_NE(p1, nullptr);
    int* p2 = pa.allocate(1);
    ASSERT_NE(p2, nullptr);
    EXPECT_NE(p2, p1 + 4);

    PageAllocator pa2(3);
    int* a = pa2.allocate(2);
    ASSERT_NE(a, nullptr);
    int* b = pa2.allocate(2);
    ASSERT_NE(b, nullptr);
    EXPECT_NE(b, a + 2);
}

TEST(PageAllocator_ResetAndPages, ClearAndReset)
{
    PageAllocator pa(4);
    pa.allocate(3);
    pa.allocate(3);
    ASSERT_EQ(pa.pages.size(), 2u);

    pa.reset();
    EXPECT_EQ(pa.pages.size(), 1u);
    EXPECT_EQ(pa.pages[0].offset, 0);

    pa.pages.clear();
    int* p = pa.allocate(1);
    ASSERT_NE(p, nullptr);
    EXPECT_EQ(pa.pages.back().offset, 1);
}

TEST(NeighborList_InitializeAndReset, Behavior)
{
    NeighborList nl;
    nl.initialize(0, 16);
    EXPECT_EQ(nl.nlocal, 0);
    EXPECT_EQ(nl.numneigh.size(), 0u);
    EXPECT_EQ(nl.firstneigh.size(), 0u);

    nl.initialize(5, 8);
    EXPECT_EQ(nl.nlocal, 5);
    EXPECT_EQ(nl.numneigh.size(), 5u);
    EXPECT_EQ(nl.firstneigh.size(), 5u);
    for (int i = 0; i < 5; ++i) {
        EXPECT_EQ(nl.numneigh[i], 0);
        EXPECT_EQ(nl.firstneigh[i], nullptr);
    }
    EXPECT_EQ(nl.allocator.pgsize, 8);

    int* storage = nl.allocator.allocate(3);
    ASSERT_NE(storage, nullptr);
    nl.numneigh[2] = 3;
    nl.firstneigh[2] = storage;

    nl.reset();
    EXPECT_EQ(nl.allocator.pages.size(), 1u);
    EXPECT_EQ(nl.allocator.pages[0].offset, 0);
    EXPECT_EQ(nl.numneigh.size(), 5u);
}


