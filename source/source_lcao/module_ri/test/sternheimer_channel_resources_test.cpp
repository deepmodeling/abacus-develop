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
    EXPECT_EQ(plan.effective_workers, 1);
}

TEST(SternheimerChannelResources, PartialOuterTeamFallsBackToNestedGridParallelism)
{
    const ModuleRI::SternheimerMemorySnapshot memory{
        ModuleRI::SternheimerMemoryAccountingMode::available,
        100000,
        0,
        1,
        "proc_meminfo"};
    const auto plan = ModuleRI::plan_sternheimer_channel_workers(32, 30, 1, 4, memory);
    EXPECT_EQ(plan.automatic_workers, 30);
    EXPECT_EQ(plan.effective_workers, 1);
}

TEST(SternheimerChannelResources, ClampsToChannelsAndOpenMPThreads)
{
    const ModuleRI::SternheimerMemorySnapshot memory{
        ModuleRI::SternheimerMemoryAccountingMode::available,
        100000,
        0,
        1,
        "proc_meminfo"};
    EXPECT_EQ(ModuleRI::plan_sternheimer_channel_workers(3, 30, 1, 0, memory).effective_workers, 1);
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

TEST(SternheimerChannelResources, FormatsWorkerDecisionDiagnostic)
{
    const ModuleRI::SternheimerMemorySnapshot memory{
        ModuleRI::SternheimerMemoryAccountingMode::node_aggregate,
        100000,
        32000,
        2,
        "cgroup_v2"};
    const auto plan = ModuleRI::plan_sternheimer_channel_workers(40, 30, 1, 4, memory);
    EXPECT_EQ(ModuleRI::format_sternheimer_channel_worker_diagnostic(memory, plan, 1, 4),
              "resource_source=cgroup_v2 accounting_mode=node_aggregate "
              "node_memory_limit_bytes=100000 memory_current_bytes=32000 "
              "memory_target_bytes=75000 local_mpi_ranks=2 grid_size=1 "
              "memory_per_worker_bytes=1920 automatic_workers=11 user_cap=4 effective_workers=1");
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

TEST(SternheimerChannelResources, ParsesCgroupSlurmAndProcMemoryUnits)
{
    using ModuleRI::detail::parse_sternheimer_kib_field;
    using ModuleRI::detail::parse_sternheimer_memory_bytes;
    using ModuleRI::detail::parse_sternheimer_slurm_mem_per_node;

    EXPECT_EQ(parse_sternheimer_memory_bytes("4096"), 4096U);
    EXPECT_EQ(parse_sternheimer_memory_bytes("110610M"), 110610ULL * 1024ULL * 1024ULL);
    EXPECT_EQ(parse_sternheimer_slurm_mem_per_node("110610"), 110610ULL * 1024ULL * 1024ULL);
    EXPECT_FALSE(parse_sternheimer_memory_bytes("max").has_value());
    EXPECT_FALSE(parse_sternheimer_memory_bytes("12MB").has_value());
    EXPECT_EQ(parse_sternheimer_kib_field("MemTotal: 4096 kB\nMemAvailable: 2048 kB\n",
                                         "MemAvailable"),
              2048ULL * 1024ULL);
    EXPECT_EQ(parse_sternheimer_kib_field("VmPeak:\t2048 kB\nVmRSS:\t1024 kB\n", "VmRSS"),
              1024ULL * 1024ULL);
    EXPECT_FALSE(parse_sternheimer_kib_field("VmRSS: invalid kB\n", "VmRSS").has_value());
}

TEST(SternheimerChannelResources, ParsesCgroupMembershipPaths)
{
    EXPECT_EQ(ModuleRI::detail::parse_sternheimer_cgroup_v2_path("0::/slurm/job_7/step_0\n"),
              "/slurm/job_7/step_0");
    EXPECT_EQ(ModuleRI::detail::parse_sternheimer_cgroup_v1_memory_path(
                  "4:cpu:/slurm/job_7\n5:memory:/slurm/job_7/step_0\n"),
              "/slurm/job_7/step_0");
    EXPECT_FALSE(ModuleRI::detail::parse_sternheimer_cgroup_v2_path("5:memory:/job\n").has_value());
}

TEST(SternheimerChannelResources, SelectsSmallestEnforcedLimitWithCgroupUsage)
{
    ModuleRI::detail::SternheimerMemoryCandidates values;
    values.cgroup_limit_bytes = 100000;
    values.cgroup_current_bytes = 20000;
    values.slurm_limit_bytes = 90000;
    values.cgroup_source = "cgroup_v2";
    const auto result = ModuleRI::detail::select_sternheimer_memory_snapshot(values, 2);
    EXPECT_EQ(result.mode, ModuleRI::SternheimerMemoryAccountingMode::node_aggregate);
    EXPECT_EQ(result.limit_bytes, 90000U);
    EXPECT_EQ(result.current_bytes, 20000U);
    EXPECT_EQ(result.local_mpi_ranks, 2);
    EXPECT_EQ(result.source, "cgroup_v2+slurm");
}

TEST(SternheimerChannelResources, UsesPerRankRssWithSlurmLimit)
{
    ModuleRI::detail::SternheimerMemoryCandidates values;
    values.slurm_limit_bytes = 100000;
    values.process_rss_bytes = 10000;
    const auto result = ModuleRI::detail::select_sternheimer_memory_snapshot(values, 2);
    EXPECT_EQ(result.mode, ModuleRI::SternheimerMemoryAccountingMode::per_rank);
    EXPECT_EQ(result.limit_bytes, 100000U);
    EXPECT_EQ(result.current_bytes, 10000U);
    EXPECT_EQ(result.source, "slurm+proc_status");
}

TEST(SternheimerChannelResources, UsesAvailableMemoryWithoutSubtractingRssAgain)
{
    ModuleRI::detail::SternheimerMemoryCandidates values;
    values.mem_available_bytes = 80000;
    values.process_rss_bytes = 10000;
    const auto result = ModuleRI::detail::select_sternheimer_memory_snapshot(values, 1);
    EXPECT_EQ(result.mode, ModuleRI::SternheimerMemoryAccountingMode::available);
    EXPECT_EQ(result.limit_bytes, 80000U);
    EXPECT_EQ(result.current_bytes, 0U);
    EXPECT_EQ(result.source, "proc_meminfo");
}

TEST(SternheimerChannelResources, TreatsV1PhysicalMemorySentinelAsUnlimited)
{
    ModuleRI::detail::SternheimerMemoryCandidates values;
    values.cgroup_limit_bytes = 0x7ffffffffffff000ULL;
    values.physical_memory_bytes = 128ULL << 30;
    values.mem_available_bytes = 64ULL << 30;
    values.cgroup_source = "cgroup_v1";
    const auto result = ModuleRI::detail::select_sternheimer_memory_snapshot(values, 1);
    EXPECT_EQ(result.mode, ModuleRI::SternheimerMemoryAccountingMode::available);
    EXPECT_EQ(result.limit_bytes, 64ULL << 30);
}

TEST(SternheimerChannelResources, FallsBackWhenNoMemorySourceIsTrustworthy)
{
    const ModuleRI::detail::SternheimerMemoryCandidates values{};
    const auto result = ModuleRI::detail::select_sternheimer_memory_snapshot(values, 3);
    EXPECT_EQ(result.mode, ModuleRI::SternheimerMemoryAccountingMode::fallback_one);
    EXPECT_EQ(result.local_mpi_ranks, 3);
    EXPECT_EQ(result.source, "unavailable");
}
