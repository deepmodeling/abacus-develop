#include "../sph_bessel_recursive.h"
#include <gtest/gtest.h>

#include <array>
#include <thread>
#include <vector>

namespace
{
typedef std::vector<std::vector<std::vector<double>>> BesselValues;

void grow_table(ModuleBase::Sph_Bessel_Recursive::D2& table)
{
    table.set_dx(0.0001);
    table.cal_jlx(7, 31, 39);
    table.cal_jlx(3, 51, 27);
    table.cal_jlx(5, 21, 61);
    table.cal_jlx(9, 33, 43);
    table.cal_jlx(2, 9, 11);
}
} // namespace

TEST(SphBesselWorkspace, EmptyAndRepeatedRelease)
{
    ModuleBase::Sph_Bessel_Recursive::D2 table;
    const std::size_t empty_bytes = table.release_workspace();
    EXPECT_EQ(empty_bytes, 0U);
    const BesselValues& empty_values = table.get_jlx();
    EXPECT_TRUE(empty_values.empty());
    grow_table(table);
    const std::size_t released_bytes = table.release_workspace();
    EXPECT_GT(released_bytes, 0U);
    const std::size_t repeated_bytes = table.release_workspace();
    EXPECT_EQ(repeated_bytes, 0U);
}

TEST(SphBesselWorkspace, ExactRaggedRebuildAndSubsequentGrowth)
{
    ModuleBase::Sph_Bessel_Recursive::D2 table;
    ModuleBase::Sph_Bessel_Recursive::D2 control;
    grow_table(table);
    grow_table(control);
    const BesselValues expected = control.get_jlx();
    for (int cycle = 0; cycle < 3; ++cycle)
    {
        const std::size_t released_bytes = table.release_workspace();
        EXPECT_GT(released_bytes, 0U);
        const BesselValues& restored = table.get_jlx();
        EXPECT_EQ(restored, expected);
    }
    table.release_workspace();
    table.cal_jlx(11, 65, 75);
    control.cal_jlx(11, 65, 75);
    const BesselValues& actual = table.get_jlx();
    const BesselValues& grown_expected = control.get_jlx();
    EXPECT_EQ(actual, grown_expected);
}

TEST(SphBesselWorkspace, ConcurrentReadAfterRelease)
{
    ModuleBase::Sph_Bessel_Recursive::D2 table;
    grow_table(table);
    const BesselValues expected = table.get_jlx();
    table.release_workspace();
    std::array<int, 4> matches = {{0, 0, 0, 0}};
    std::vector<std::thread> readers;
    for (int index = 0; index < 4; ++index)
    {
        readers.emplace_back([&table, &expected, &matches, index]() {
            const BesselValues& actual = table.get_jlx();
            matches[index] = actual == expected;
        });
    }
    for (std::thread& reader : readers)
    {
        reader.join();
    }
    for (const int match : matches)
    {
        EXPECT_EQ(match, 1);
    }
}

TEST(SphBesselWorkspace, CopyAndMoveAfterRelease)
{
    std::vector<ModuleBase::Sph_Bessel_Recursive::D2> tables(1);
    grow_table(tables[0]);
    const BesselValues expected = tables[0].get_jlx();
    tables[0].release_workspace();
    ModuleBase::Sph_Bessel_Recursive::D2 copied = tables[0];
    tables.emplace_back();
    const BesselValues& moved_values = tables[0].get_jlx();
    const BesselValues& copied_values = copied.get_jlx();
    EXPECT_EQ(moved_values, expected);
    EXPECT_EQ(copied_values, expected);
}
