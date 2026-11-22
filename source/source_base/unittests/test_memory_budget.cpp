#include "source_base/memory_budget.h"

#include <gtest/gtest.h>

#include <cstddef>
#include <fstream>
#include <limits>
#include <string>
#include <sys/stat.h>

namespace
{

void write_file(const std::string& path, const std::string& content)
{
    std::ofstream ofs(path);
    ASSERT_TRUE(ofs);
    ofs << content;
}

std::string make_v2_fixture()
{
    const std::string dir = testing::TempDir() + "/abacus_memory_budget_v2";
    mkdir(dir.c_str(), 0700);
    write_file(dir + "/memory.max", "1536000000000\n");
    write_file(dir + "/memory.current", "100000000000\n");
    write_file(dir + "/memory.stat", "inactive_file 20000000000\n");
    return dir;
}

std::string make_v1_fixture()
{
    const std::string dir = testing::TempDir() + "/abacus_memory_budget_v1";
    mkdir(dir.c_str(), 0700);
    mkdir((dir + "/memory").c_str(), 0700);
    write_file(dir + "/memory/memory.limit_in_bytes", "1536000000000\n");
    write_file(dir + "/memory/memory.usage_in_bytes", "100000000000\n");
    write_file(dir + "/memory/memory.stat", "total_inactive_file 20000000000\n");
    return dir;
}

} // namespace

TEST(MemoryBudget, CgroupV2ReadsLimitAndCurrent)
{
    const ModuleBase::MemoryStatus status = ModuleBase::get_memory_status(make_v2_fixture());

    EXPECT_TRUE(status.from_cgroup);
    EXPECT_EQ(status.limit_bytes, 1536000000000ULL);
    EXPECT_EQ(status.current_bytes, 100000000000ULL);
    EXPECT_EQ(status.reclaimable_bytes, 20000000000ULL);
    EXPECT_EQ(status.available_bytes, 1456000000000ULL);
}

TEST(MemoryBudget, CgroupV1ReadsLimitAndCurrent)
{
    const ModuleBase::MemoryStatus status = ModuleBase::get_memory_status(make_v1_fixture());

    EXPECT_TRUE(status.from_cgroup);
    EXPECT_EQ(status.limit_bytes, 1536000000000ULL);
    EXPECT_EQ(status.current_bytes, 100000000000ULL);
    EXPECT_EQ(status.reclaimable_bytes, 20000000000ULL);
    EXPECT_EQ(status.available_bytes, 1456000000000ULL);
}

TEST(MemoryBudget, CgroupUnlimitedKeepsCurrentButNoLimit)
{
    const std::string dir = make_v2_fixture();
    write_file(dir + "/memory.max", "max\n");

    const ModuleBase::MemoryStatus status = ModuleBase::get_memory_status(dir);
    EXPECT_TRUE(status.from_cgroup);
    EXPECT_EQ(status.limit_bytes, 0ULL);
    EXPECT_EQ(status.current_bytes, 100000000000ULL);
}

TEST(MemoryBudget, EstimatorMatchesIssue8342Dimensions)
{
    // NLOCAL ~ 377884 split over 56 ranks gives a local block of ~2.55e9
    // elements; gamma-only genelpa keeps five such dense matrices.
    const std::size_t nrloc = 50500;
    const std::size_t ncloc = 50500;
    const std::size_t local = nrloc * ncloc;
    const ModuleBase::LcaoDenseEstimate estimate
        = ModuleBase::estimate_lcao_dense(nrloc, ncloc, 73820, 1, 1, true, true, true, 8);

    EXPECT_EQ(estimate.hk_bytes, local * 8);
    EXPECT_EQ(estimate.sk_bytes, local * 8);
    EXPECT_EQ(estimate.dmk_bytes, local * 8);
    EXPECT_EQ(estimate.psi_bytes, local * 8);
    EXPECT_EQ(estimate.elpa_work_bytes, local * 8);
    EXPECT_EQ(estimate.per_rank_total_bytes, local * 8 * 5);
    EXPECT_GT(estimate.per_rank_total_bytes * 56, 1500000000000ULL);
}

TEST(MemoryBudget, EstimatorHandlesMultiKAndNonGamma)
{
    const ModuleBase::LcaoDenseEstimate estimate
        = ModuleBase::estimate_lcao_dense(1000, 1000, 100, 4, 1, false, false, false, 8);

    EXPECT_EQ(estimate.hk_bytes, 8 * 1000 * 1000);
    EXPECT_EQ(estimate.sk_bytes, 8 * 1000 * 1000);
    EXPECT_EQ(estimate.dmk_bytes, 4 * 8 * 1000 * 1000);
    EXPECT_EQ(estimate.psi_bytes, 4 * 100 * 1000 * 8);
    EXPECT_EQ(estimate.elpa_work_bytes, 0ULL);
}

TEST(MemoryBudget, EstimatorSaturatesOnOverflow)
{
    const ModuleBase::LcaoDenseEstimate estimate
        = ModuleBase::estimate_lcao_dense(std::numeric_limits<std::size_t>::max(),
                                          std::numeric_limits<std::size_t>::max(),
                                          1,
                                          1,
                                          1,
                                          true,
                                          true,
                                          false,
                                          16);
    EXPECT_EQ(estimate.hk_bytes, std::numeric_limits<std::size_t>::max());
    EXPECT_EQ(estimate.per_rank_total_bytes, std::numeric_limits<std::size_t>::max());
}

TEST(MemoryBudget, DecideAbortsWhenLowerBoundDoesNotFit)
{
    EXPECT_EQ(ModuleBase::decide_memory(1600000000000ULL, 1800000000000ULL, 1500000000000ULL),
              ModuleBase::MemoryDecision::Abort);
}

TEST(MemoryBudget, DecideWarnsWhenProjectionDoesNotFit)
{
    EXPECT_EQ(ModuleBase::decide_memory(1000000000000ULL, 1600000000000ULL, 1500000000000ULL),
              ModuleBase::MemoryDecision::Warn);
}

TEST(MemoryBudget, DecideOkWhenProjectionFits)
{
    EXPECT_EQ(ModuleBase::decide_memory(500000000000ULL, 600000000000ULL, 1500000000000ULL),
              ModuleBase::MemoryDecision::Ok);
}
