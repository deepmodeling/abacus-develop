#include "source_lcao/module_ri/sternheimer_channel_resources.h"

#include <complex>
#include <gtest/gtest.h>
#include <limits>
#include <stdexcept>

TEST(SternheimerChannelResources, EstimatesOneHundredTwentyComplexGridVectors)
{
    EXPECT_EQ(ModuleRI::estimate_sternheimer_channel_worker_bytes(1000),
              120ULL * 1000ULL * sizeof(std::complex<double>));
}

TEST(SternheimerChannelResources, PlansFromNodeAggregateMemory)
{
    const ModuleRI::SternheimerMemorySnapshot memory{
        ModuleRI::SternheimerMemoryAccountingMode::node_aggregate,
        100000,
        32000,
        2,
        "cgroup_v2"};
    const auto plan = ModuleRI::plan_sternheimer_channel_workers(40, 30, 1, 0, memory);
    EXPECT_EQ(plan.automatic_workers, 11);
    EXPECT_EQ(plan.effective_workers, 11);
}

TEST(SternheimerChannelResources, UserCapOnlyReducesAutomaticCount)
{
    const ModuleRI::SternheimerMemorySnapshot memory{
        ModuleRI::SternheimerMemoryAccountingMode::available,
        100000,
        0,
        1,
        "proc_meminfo"};
    const auto plan = ModuleRI::plan_sternheimer_channel_workers(32, 30, 1, 4, memory);
    EXPECT_EQ(plan.automatic_workers, 30);
    EXPECT_EQ(plan.effective_workers, 4);
}

TEST(SternheimerChannelResources, ClampsToChannelsAndOpenMPThreads)
{
    const ModuleRI::SternheimerMemorySnapshot memory{
        ModuleRI::SternheimerMemoryAccountingMode::available,
        100000,
        0,
        1,
        "proc_meminfo"};
    EXPECT_EQ(ModuleRI::plan_sternheimer_channel_workers(3, 30, 1, 0, memory).effective_workers, 3);
    EXPECT_EQ(ModuleRI::plan_sternheimer_channel_workers(40, 5, 1, 0, memory).effective_workers, 5);
}

TEST(SternheimerChannelResources, DividesPerRankLimitBeforeSubtractingRss)
{
    const ModuleRI::SternheimerMemorySnapshot memory{
        ModuleRI::SternheimerMemoryAccountingMode::per_rank,
        100000,
        10000,
        2,
        "slurm+proc_status"};
    const auto plan = ModuleRI::plan_sternheimer_channel_workers(40, 30, 1, 0, memory);
    EXPECT_EQ(plan.automatic_workers, 14);
}

TEST(SternheimerChannelResources, DividesAvailableMemoryAcrossLocalRanks)
{
    const ModuleRI::SternheimerMemorySnapshot memory{
        ModuleRI::SternheimerMemoryAccountingMode::available,
        100000,
        0,
        2,
        "proc_meminfo"};
    const auto plan = ModuleRI::plan_sternheimer_channel_workers(40, 30, 1, 0, memory);
    EXPECT_EQ(plan.automatic_workers, 19);
}

TEST(SternheimerChannelResources, ZeroCapUsesAutomaticCount)
{
    const ModuleRI::SternheimerMemorySnapshot memory{
        ModuleRI::SternheimerMemoryAccountingMode::available,
        100000,
        0,
        1,
        "proc_meminfo"};
    const auto plan = ModuleRI::plan_sternheimer_channel_workers(32, 30, 1, 0, memory);
    EXPECT_EQ(plan.effective_workers, plan.automatic_workers);
}

TEST(SternheimerChannelResources, RejectsDetectedBudgetBelowOneWorker)
{
    const ModuleRI::SternheimerMemorySnapshot memory{
        ModuleRI::SternheimerMemoryAccountingMode::node_aggregate,
        1000,
        0,
        1,
        "cgroup_v2"};
    EXPECT_THROW(ModuleRI::plan_sternheimer_channel_workers(4, 4, 1, 0, memory), std::runtime_error);
}

TEST(SternheimerChannelResources, UnavailableResourcesFallBackToOneWorker)
{
    const ModuleRI::SternheimerMemorySnapshot memory{};
    EXPECT_EQ(ModuleRI::plan_sternheimer_channel_workers(32, 30, 1, 0, memory).effective_workers, 1);
    EXPECT_EQ(ModuleRI::plan_sternheimer_channel_workers(32, 30, 1, 8, memory).effective_workers, 1);
}

TEST(SternheimerChannelResources, RejectsInvalidInputsAndOverflow)
{
    const ModuleRI::SternheimerMemorySnapshot memory{};
    EXPECT_THROW(ModuleRI::plan_sternheimer_channel_workers(-1, 1, 1, 0, memory), std::invalid_argument);
    EXPECT_THROW(ModuleRI::plan_sternheimer_channel_workers(1, 0, 1, 0, memory), std::invalid_argument);
    EXPECT_THROW(ModuleRI::plan_sternheimer_channel_workers(1, 1, 1, -1, memory), std::invalid_argument);
    EXPECT_THROW(ModuleRI::plan_sternheimer_channel_workers(1, 1, 0, 0, memory), std::invalid_argument);
    EXPECT_THROW(ModuleRI::estimate_sternheimer_channel_worker_bytes(std::numeric_limits<std::size_t>::max()),
                 std::overflow_error);
}
