#include "source_io/module_output/provenance_header.h"

#include "source_main/version.h"

#include "gtest/gtest.h"

#include <string>

/**
 * Tests for the provenance header shared by the band, wavefunction and
 * band-decomposed energy writers.
 */

TEST(ProvenanceHeaderTest, TimestampFormat)
{
    const std::string stamp = ModuleIO::local_timestamp();
    ASSERT_EQ(stamp.size(), 19u);
    EXPECT_EQ(stamp[4], '-');
    EXPECT_EQ(stamp[7], '-');
    EXPECT_EQ(stamp[10], ' ');
    EXPECT_EQ(stamp[13], ':');
    EXPECT_EQ(stamp[16], ':');
}

TEST(ProvenanceHeaderTest, HeaderCarriesVersionAndTimestamp)
{
    const std::string header = ModuleIO::provenance_header("2026-01-01 00:00:00");
    const std::string expected = "# ABACUS version: " + std::string(VERSION)
                                 + "\n# Written at 2026-01-01 00:00:00\n";
    EXPECT_EQ(header, expected);
}
