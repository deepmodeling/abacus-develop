#include <gtest/gtest.h>
#include "../neighbor_list.h"

TEST(PageAllocator_Constructors, DefaultAndCustom)
{
    PageAllocator pa_def;
    EXPECT_EQ(pa_def.get_pgsize(), PageAllocator::default_pgsize);

    PageAllocator pa(4);
    EXPECT_EQ(pa.get_pgsize(), 4);
}

TEST(PageAllocator_AllocateEdgeCases, ZeroNegative)
{
    PageAllocator pa(8);
    EXPECT_EQ(pa.allocate(0), nullptr);
    EXPECT_EQ(pa.allocate(-5), nullptr);
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

TEST(PageAllocator_Reset, ClearAndReset)
{
    PageAllocator pa(4);
    pa.allocate(3);
    pa.allocate(3);

    pa.reset();
    int* p = pa.allocate(1);
    ASSERT_NE(p, nullptr);
}

TEST(NeighborList_InitializeAndReset, Behavior)
{
    NeighborList nl;
    nl.initialize(0, 16);
    EXPECT_EQ(nl.get_nlocal(), 0);

    nl.initialize(5, 8);
    EXPECT_EQ(nl.get_nlocal(), 5);
    for (int i = 0; i < 5; ++i) {
        EXPECT_EQ(nl.get_numneigh(i), 0);
        EXPECT_EQ(nl.get_firstneigh(i), nullptr);
    }

    nl.reset();
    EXPECT_EQ(nl.get_nlocal(), 5);
}

TEST(NeighborList_Getters, Accessors)
{
    NeighborList nl;
    nl.initialize(3, 16);
    EXPECT_EQ(nl.get_nlocal(), 3);
    EXPECT_EQ(nl.get_numneigh(0), 0);
    EXPECT_EQ(nl.get_firstneigh(0), nullptr);
}