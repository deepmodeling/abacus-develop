#include "source_lcao/module_ri/sternheimer_channel_resources.h"
#include "source_lcao/module_ri/sternheimer_abacus_st_smoke.h"

#include <complex>
#include <gtest/gtest.h>
#include <limits>
#include <stdexcept>

TEST(SternheimerChannelResources, EstimatesOneHundredTwentyComplexGridVectors)
{
    EXPECT_EQ(ModuleRI::estimate_sternheimer_channel_worker_bytes(1000),
              120ULL * 1000ULL * sizeof(std::complex<double>));
}

TEST(SternheimerChannelResources, GroupsChannelsIntoStableMicroBatches)
{
    const auto batches = ModuleRI::make_sternheimer_channel_batches(10, 4);
    ASSERT_EQ(batches.size(), 3U);
    EXPECT_EQ(batches[0].begin, 0);
    EXPECT_EQ(batches[0].size, 4);
    EXPECT_EQ(batches[1].begin, 4);
    EXPECT_EQ(batches[1].size, 4);
    EXPECT_EQ(batches[2].begin, 8);
    EXPECT_EQ(batches[2].size, 2);

    const auto scalar = ModuleRI::make_sternheimer_channel_batches(3, 1);
    ASSERT_EQ(scalar.size(), 3U);
    EXPECT_EQ(scalar[0].begin, 0);
    EXPECT_EQ(scalar[1].begin, 1);
    EXPECT_EQ(scalar[2].begin, 2);
    EXPECT_THROW(ModuleRI::make_sternheimer_channel_batches(-1, 4), std::invalid_argument);
    EXPECT_THROW(ModuleRI::make_sternheimer_channel_batches(4, 0), std::invalid_argument);
}

TEST(SternheimerChannelResources, OwnedChannelBatchesKeepRequestedWidthOwnership)
{
    const auto first = ModuleRI::make_sternheimer_owned_channel_batches(10, 4, 2, 2, 0, 0);
    ASSERT_EQ(first.size(), 3U);
    EXPECT_EQ(first[0].begin, 0);
    EXPECT_EQ(first[0].size, 2);
    EXPECT_EQ(first[1].begin, 2);
    EXPECT_EQ(first[1].size, 2);
    EXPECT_EQ(first[2].begin, 8);
    EXPECT_EQ(first[2].size, 2);
    const auto second = ModuleRI::make_sternheimer_owned_channel_batches(10, 4, 3, 2, 1, 0);
    ASSERT_EQ(second.size(), 2U);
    EXPECT_EQ(second[0].begin, 4);
    EXPECT_EQ(second[0].size, 3);
    EXPECT_EQ(second[1].begin, 7);
    EXPECT_EQ(second[1].size, 1);
}

TEST(SternheimerChannelResources, OwnedChannelBatchesMixedWidthsCoverEveryChannelExactlyOnce)
{
    for (const int num_channels: {1, 2, 3, 7, 10, 17, 31})
    {
        for (const int ownership_width: {1, 2, 4, 7})
        {
            for (const int replicas: {1, 2, 3, 5})
            {
                for (const int occupied_band: {0, 1, 3, std::numeric_limits<int>::max()})
                {
                    SCOPED_TRACE(std::to_string(num_channels) + ":" + std::to_string(ownership_width)
                                 + ":" + std::to_string(replicas) + ":" + std::to_string(occupied_band));
                    const int ownership_batch_count = 1 + (num_channels - 1) / ownership_width;
                    std::vector<int> visits(static_cast<std::size_t>(num_channels), 0);
                    for (int replica = 0; replica != replicas; ++replica)
                    {
                        const int worker_width = 1 + replica % 3;
                        const auto batches = ModuleRI::make_sternheimer_owned_channel_batches(
                            num_channels, ownership_width, worker_width, replicas, replica, occupied_band);
                        int previous_end = 0;
                        for (const auto& batch: batches)
                        {
                            ASSERT_GE(batch.begin, previous_end);
                            ASSERT_LT(batch.begin, num_channels);
                            ASSERT_GT(batch.size, 0);
                            ASSERT_LE(batch.size, worker_width);
                            ASSERT_LE(batch.size, num_channels - batch.begin);
                            EXPECT_EQ(batch.begin / ownership_width,
                                      (batch.begin + batch.size - 1) / ownership_width);
                            for (int channel = batch.begin; channel != batch.begin + batch.size; ++channel)
                            {
                                const int owner = ModuleRI::sternheimer_channel_batch_replica_owner(
                                    occupied_band, channel / ownership_width, ownership_batch_count, replicas);
                                EXPECT_EQ(owner, replica);
                                ++visits[static_cast<std::size_t>(channel)];
                            }
                            previous_end = batch.begin + batch.size;
                        }
                    }
                    for (const int count: visits)
                    {
                        EXPECT_EQ(count, 1);
                    }
                }
            }
        }
    }
}

TEST(SternheimerChannelResources, OwnedChannelBatchesEqualWidthsPreserveLegacyBandRotation)
{
    for (const int width: {1, 2, 4})
    {
        const auto global_batches = ModuleRI::make_sternheimer_channel_batches(11, width);
        for (int occupied_band = 0; occupied_band != 7; ++occupied_band)
        {
            for (int replica = 0; replica != 3; ++replica)
            {
                std::vector<ModuleRI::SternheimerChannelBatch> expected;
                for (std::size_t chunk = 0; chunk != global_batches.size(); ++chunk)
                {
                    if (ModuleRI::sternheimer_channel_batch_replica_owner(
                            occupied_band, static_cast<int>(chunk), static_cast<int>(global_batches.size()), 3)
                        == replica)
                    {
                        expected.push_back(global_batches[chunk]);
                    }
                }
                const auto actual = ModuleRI::make_sternheimer_owned_channel_batches(
                    11, width, width, 3, replica, occupied_band);
                ASSERT_EQ(actual.size(), expected.size());
                for (std::size_t batch = 0; batch != expected.size(); ++batch)
                {
                    EXPECT_EQ(actual[batch].begin, expected[batch].begin);
                    EXPECT_EQ(actual[batch].size, expected[batch].size);
                }
            }
        }
    }
}

TEST(SternheimerChannelResources, OwnedChannelBatchesAllowEmptyChannelsAndIdleReplicas)
{
    EXPECT_TRUE(ModuleRI::make_sternheimer_owned_channel_batches(0, 4, 2, 3, 1, 0).empty());
    EXPECT_TRUE(ModuleRI::make_sternheimer_owned_channel_batches(3, 4, 2, 8, 7, 0).empty());
    EXPECT_TRUE(ModuleRI::make_sternheimer_owned_channel_batches(3, 4, 2, 8, 0, 7).empty());
    const auto rotated = ModuleRI::make_sternheimer_owned_channel_batches(3, 4, 8, 8, 7, 7);
    ASSERT_EQ(rotated.size(), 1U);
    EXPECT_EQ(rotated[0].begin, 0);
    EXPECT_EQ(rotated[0].size, 3);
}

TEST(SternheimerChannelResources, OwnedChannelBatchesRejectInvalidInputs)
{
    EXPECT_THROW(ModuleRI::make_sternheimer_owned_channel_batches(-1, 2, 1, 2, 0, 0), std::invalid_argument);
    for (const int invalid: {0, -1})
    {
        EXPECT_THROW(ModuleRI::make_sternheimer_owned_channel_batches(4, invalid, 1, 2, 0, 0), std::invalid_argument);
        EXPECT_THROW(ModuleRI::make_sternheimer_owned_channel_batches(4, 2, invalid, 2, 0, 0), std::invalid_argument);
        EXPECT_THROW(ModuleRI::make_sternheimer_owned_channel_batches(4, 2, 1, invalid, 0, 0), std::invalid_argument);
    }
    for (const int invalid: {-1, 2})
    {
        EXPECT_THROW(ModuleRI::make_sternheimer_owned_channel_batches(4, 2, 1, 2, invalid, 0), std::invalid_argument);
    }
    EXPECT_THROW(ModuleRI::make_sternheimer_owned_channel_batches(4, 2, 1, 2, 0, -1), std::invalid_argument);
    EXPECT_THROW(ModuleRI::make_sternheimer_owned_channel_batches(0, 0, 1, 2, 0, 0), std::invalid_argument);
}

TEST(SternheimerChannelResources, OwnedChannelBatchesAvoidIntegerOverflow)
{
    constexpr int maximum = std::numeric_limits<int>::max();
    const auto wide = ModuleRI::make_sternheimer_owned_channel_batches(maximum, maximum, maximum,
                                                                      maximum, maximum - 1, maximum - 1);
    ASSERT_EQ(wide.size(), 1U);
    EXPECT_EQ(wide[0].begin, 0);
    EXPECT_EQ(wide[0].size, maximum);
    const auto tail = ModuleRI::make_sternheimer_owned_channel_batches(maximum, maximum - 1, maximum, 2, 1, 0);
    ASSERT_EQ(tail.size(), 1U);
    EXPECT_EQ(tail[0].begin, maximum - 1);
    EXPECT_EQ(tail[0].size, 1);
    const auto sparse = ModuleRI::make_sternheimer_owned_channel_batches(maximum, 1, 1,
                                                                        maximum, maximum - 1, maximum);
    ASSERT_EQ(sparse.size(), 1U);
    EXPECT_EQ(sparse[0].begin, maximum - 1);
    EXPECT_EQ(sparse[0].size, 1);
}

TEST(SternheimerChannelResources, BatchWorkerPlanUsesEveryAvailableOuterThread)
{
    const ModuleRI::SternheimerMemorySnapshot memory{ModuleRI::SternheimerMemoryAccountingMode::available,
                                                     100000000,
                                                     0,
                                                     1,
                                                     "test"};
    const auto scalar = ModuleRI::plan_sternheimer_owned_channel_workers(100, 100, 30, 100, 0, memory, 1);
    const auto batch = ModuleRI::plan_sternheimer_owned_channel_workers(100, 100, 30, 100, 0, memory, 4);

    EXPECT_EQ(batch.channel_batch_width, 4);
    EXPECT_EQ(batch.batch_tasks, 25);
    EXPECT_EQ(batch.memory_per_worker_bytes, 4 * scalar.memory_per_worker_bytes);
    EXPECT_EQ(scalar.automatic_workers, 30);
    EXPECT_EQ(batch.automatic_workers, 25);
    EXPECT_EQ(batch.effective_workers, 25);
    EXPECT_LE(static_cast<std::uint64_t>(batch.automatic_workers) * batch.memory_per_worker_bytes,
              batch.increment_bytes_per_rank);
}

TEST(SternheimerChannelResources, PlansFromNodeAggregateMemory)
{
    const ModuleRI::SternheimerMemorySnapshot memory{ModuleRI::SternheimerMemoryAccountingMode::node_aggregate,
                                                     100000,
                                                     32000,
                                                     2,
                                                     "cgroup_v2"};
    const auto plan = ModuleRI::plan_sternheimer_channel_workers(40, 30, 1, 0, memory);
    EXPECT_EQ(plan.automatic_workers, 11);
    EXPECT_EQ(plan.effective_workers, 11);
}

TEST(SternheimerChannelResources, PartialOuterTeamFallsBackToNestedGridParallelism)
{
    const ModuleRI::SternheimerMemorySnapshot memory{ModuleRI::SternheimerMemoryAccountingMode::available,
                                                     100000,
                                                     0,
                                                     1,
                                                     "proc_meminfo"};
    const auto plan = ModuleRI::plan_sternheimer_channel_workers(32, 30, 1, 4, memory);
    EXPECT_EQ(plan.automatic_workers, 30);
    EXPECT_EQ(plan.effective_workers, 1);
}

TEST(SternheimerChannelResources, BatchColumnsDoNotCountAsOuterWorkerThreads)
{
    const ModuleRI::SternheimerMemorySnapshot memory{ModuleRI::SternheimerMemoryAccountingMode::available,
                                                     100000000,
                                                     0,
                                                     1,
                                                     "proc_meminfo"};
    const auto plan = ModuleRI::plan_sternheimer_channel_workers(8, 30, 100, 0, memory, 4);
    EXPECT_EQ(plan.automatic_workers, 2);
    EXPECT_EQ(plan.effective_workers, 1);
}

TEST(SternheimerChannelResources, DownshiftsRequestedBatchWidthToDetectedMemoryBudget)
{
    const ModuleRI::SternheimerMemorySnapshot memory{ModuleRI::SternheimerMemoryAccountingMode::node_aggregate,
                                                     100000,
                                                     50000,
                                                     1,
                                                     "cgroup_v1"};
    const auto plan = ModuleRI::plan_sternheimer_channel_workers(40, 30, 1, 0, memory, 16);

    EXPECT_EQ(plan.channel_batch_width, 13);
    EXPECT_EQ(plan.batch_tasks, 4);
    EXPECT_EQ(plan.memory_per_worker_bytes, 13 * 120 * sizeof(std::complex<double>));
    EXPECT_LE(plan.memory_per_worker_bytes, plan.increment_bytes_per_rank);
    EXPECT_EQ(plan.automatic_workers, 1);
    EXPECT_EQ(plan.effective_workers, 1);
}

TEST(SternheimerChannelResources, KeepsOuterParallelismWhenMostThreadsCanWork)
{
    const ModuleRI::SternheimerMemorySnapshot memory{ModuleRI::SternheimerMemoryAccountingMode::available,
                                                     56320,
                                                     0,
                                                     1,
                                                     "proc_meminfo"};
    const auto plan = ModuleRI::plan_sternheimer_channel_workers(40, 30, 1, 0, memory);
    EXPECT_EQ(plan.automatic_workers, 22);
    EXPECT_EQ(plan.effective_workers, 22);
}

TEST(SternheimerChannelResources, ClampsToChannelsAndOpenMPThreads)
{
    const ModuleRI::SternheimerMemorySnapshot memory{ModuleRI::SternheimerMemoryAccountingMode::available,
                                                     100000,
                                                     0,
                                                     1,
                                                     "proc_meminfo"};
    EXPECT_EQ(ModuleRI::plan_sternheimer_channel_workers(3, 30, 1, 0, memory).effective_workers, 1);
    EXPECT_EQ(ModuleRI::plan_sternheimer_channel_workers(40, 5, 1, 0, memory).effective_workers, 5);
}

TEST(SternheimerChannelResources, DividesPerRankLimitBeforeSubtractingRss)
{
    const ModuleRI::SternheimerMemorySnapshot memory{ModuleRI::SternheimerMemoryAccountingMode::per_rank,
                                                     100000,
                                                     10000,
                                                     2,
                                                     "slurm+proc_status"};
    const auto plan = ModuleRI::plan_sternheimer_channel_workers(40, 30, 1, 0, memory);
    EXPECT_EQ(plan.automatic_workers, 14);
}

TEST(SternheimerChannelResources, DividesAvailableMemoryAcrossLocalRanks)
{
    const ModuleRI::SternheimerMemorySnapshot memory{ModuleRI::SternheimerMemoryAccountingMode::available,
                                                     100000,
                                                     0,
                                                     2,
                                                     "proc_meminfo"};
    const auto plan = ModuleRI::plan_sternheimer_channel_workers(40, 30, 1, 0, memory);
    EXPECT_EQ(plan.automatic_workers, 19);
}

TEST(SternheimerChannelResources, ZeroCapUsesAutomaticCount)
{
    const ModuleRI::SternheimerMemorySnapshot memory{ModuleRI::SternheimerMemoryAccountingMode::available,
                                                     100000,
                                                     0,
                                                     1,
                                                     "proc_meminfo"};
    const auto plan = ModuleRI::plan_sternheimer_channel_workers(32, 30, 1, 0, memory);
    EXPECT_EQ(plan.effective_workers, plan.automatic_workers);
}

TEST(SternheimerChannelResources, FormatsWorkerDecisionDiagnostic)
{
    const ModuleRI::SternheimerMemorySnapshot memory{ModuleRI::SternheimerMemoryAccountingMode::node_aggregate,
                                                     100000,
                                                     32000,
                                                     2,
                                                     "cgroup_v2"};
    const auto plan = ModuleRI::plan_sternheimer_channel_workers(40, 30, 1, 4, memory);
    EXPECT_EQ(ModuleRI::format_sternheimer_channel_worker_diagnostic(memory, plan, 1, 4),
              "resource_source=cgroup_v2 accounting_mode=node_aggregate "
              "node_memory_limit_bytes=100000 memory_current_bytes=32000 "
              "memory_target_bytes=75000 local_mpi_ranks=2 grid_size=1 "
              "memory_per_worker_bytes=1920 channel_batch_width=1 batch_tasks=40 "
              "automatic_workers=11 user_cap=4 effective_workers=1");
}

TEST(SternheimerChannelResources, RejectsDetectedBudgetBelowOneWorker)
{
    const ModuleRI::SternheimerMemorySnapshot memory{ModuleRI::SternheimerMemoryAccountingMode::node_aggregate,
                                                     1000,
                                                     0,
                                                     1,
                                                     "cgroup_v2"};
    EXPECT_THROW(ModuleRI::plan_sternheimer_channel_workers(4, 4, 1, 0, memory), std::runtime_error);
}

TEST(SternheimerChannelResources, BaselineRejectionReportsSnapshotAndPreservesSeventyFivePercentTarget)
{
    for (const auto mode: {ModuleRI::SternheimerMemoryAccountingMode::node_aggregate,
                          ModuleRI::SternheimerMemoryAccountingMode::per_rank})
    {
        const std::uint64_t target
            = mode == ModuleRI::SternheimerMemoryAccountingMode::node_aggregate ? 75000 : 37500;
        for (const std::uint64_t current: {target, target + 1})
        {
            SCOPED_TRACE(std::to_string(current));
            const ModuleRI::SternheimerMemorySnapshot memory{mode, 100000, current, 2, "test_limit"};
            try
            {
                (void)ModuleRI::plan_sternheimer_channel_workers(4, 4, 1, 1, memory);
                FAIL() << "A one-worker cap must not bypass the 75 percent target.";
            }
            catch (const std::runtime_error& error)
            {
                const std::string message = error.what();
                EXPECT_NE(message.find("75 percent target"), std::string::npos);
                EXPECT_NE(message.find("memory_current_bytes=" + std::to_string(current)), std::string::npos);
                EXPECT_NE(message.find("memory_target_bytes=" + std::to_string(target)), std::string::npos);
                EXPECT_NE(message.find("node_memory_limit_bytes=100000"), std::string::npos);
                EXPECT_NE(message.find("resource_source=test_limit"), std::string::npos);
                EXPECT_NE(message.find("accounting_mode=" + ModuleRI::sternheimer_memory_accounting_mode_name(mode)),
                          std::string::npos);
                EXPECT_NE(message.find("local_mpi_ranks=2"), std::string::npos);
            }
        }
    }
}

TEST(SternheimerChannelResources, InsufficientWorkerHeadroomReportsCurrentTargetAndLimit)
{
    const ModuleRI::SternheimerMemorySnapshot memory{ModuleRI::SternheimerMemoryAccountingMode::node_aggregate,
                                                     100000,
                                                     74999,
                                                     1,
                                                     "cgroup_v2"};
    try
    {
        (void)ModuleRI::plan_sternheimer_channel_workers(4, 4, 1, 0, memory);
        FAIL() << "One byte of headroom must not admit a worker.";
    }
    catch (const std::runtime_error& error)
    {
        const std::string message = error.what();
        EXPECT_NE(message.find("75 percent target"), std::string::npos);
        EXPECT_NE(message.find("memory_current_bytes=74999"), std::string::npos);
        EXPECT_NE(message.find("memory_target_bytes=75000"), std::string::npos);
        EXPECT_NE(message.find("node_memory_limit_bytes=100000"), std::string::npos);
    }
}

TEST(SternheimerChannelResources, ExactWorkerHeadroomRemainsAdmissible)
{
    const auto worker_bytes = ModuleRI::estimate_sternheimer_channel_worker_bytes(1);
    ModuleRI::SternheimerMemorySnapshot memory{ModuleRI::SternheimerMemoryAccountingMode::node_aggregate,
                                               100000,
                                               75000 - worker_bytes,
                                               1,
                                               "cgroup_v2"};
    const auto plan = ModuleRI::plan_sternheimer_channel_workers(4, 4, 1, 1, memory);
    EXPECT_EQ(plan.target_bytes, 75000U);
    EXPECT_EQ(plan.increment_bytes_per_rank, worker_bytes);
    EXPECT_EQ(plan.effective_workers, 1);
    ++memory.current_bytes;
    EXPECT_THROW(ModuleRI::plan_sternheimer_channel_workers(4, 4, 1, 1, memory), std::runtime_error);
}

TEST(SternheimerChannelResources, BestEffortHeapTrimPreservesLiveDataAndDoesNotOverrideBudget)
{
    static_assert(noexcept(ModuleRI::trim_sternheimer_process_heap()), "Heap trimming must be best effort.");
    const std::complex<double> value(3.25, -0.75);
    const std::vector<std::complex<double>> live_values(1024, value);
    const ModuleRI::SternheimerMemorySnapshot memory{ModuleRI::SternheimerMemoryAccountingMode::node_aggregate,
                                                     100000,
                                                     75000,
                                                     1,
                                                     "cgroup_v2"};
    EXPECT_THROW(ModuleRI::plan_sternheimer_channel_workers(4, 4, 1, 1, memory), std::runtime_error);
    for (int attempt = 0; attempt != 2; ++attempt)
    {
        const auto status = ModuleRI::trim_sternheimer_process_heap();
#if defined(__GLIBC__)
        EXPECT_TRUE(status == ModuleRI::SternheimerHeapTrimStatus::no_pages_released
                    || status == ModuleRI::SternheimerHeapTrimStatus::pages_released);
#else
        EXPECT_EQ(status, ModuleRI::SternheimerHeapTrimStatus::unsupported);
#endif
        for (const auto& element: live_values)
        {
            EXPECT_EQ(element, value);
        }
        EXPECT_THROW(ModuleRI::plan_sternheimer_channel_workers(4, 4, 1, 1, memory), std::runtime_error);
    }
}

TEST(SternheimerChannelResources, ReservesSharedMemoryAcrossNodeAggregateRanks)
{
    const ModuleRI::SternheimerMemorySnapshot memory{ModuleRI::SternheimerMemoryAccountingMode::node_aggregate,
                                                     100000,
                                                     32000,
                                                     2,
                                                     "cgroup_v2"};
    const auto adjusted = ModuleRI::reserve_sternheimer_shared_memory(memory, 10000);
    EXPECT_EQ(adjusted.current_bytes, 52000U);
    EXPECT_EQ(adjusted.limit_bytes, memory.limit_bytes);
    EXPECT_EQ(adjusted.mode, memory.mode);
    EXPECT_EQ(adjusted.local_mpi_ranks, memory.local_mpi_ranks);
    EXPECT_EQ(adjusted.source, memory.source);
    EXPECT_EQ(memory.current_bytes, 32000U);
    EXPECT_EQ(ModuleRI::plan_sternheimer_channel_workers(40, 30, 1, 0, adjusted).automatic_workers, 5);
}

TEST(SternheimerChannelResources, ReservesSharedMemoryOnceForPerRankAccounting)
{
    const ModuleRI::SternheimerMemorySnapshot memory{ModuleRI::SternheimerMemoryAccountingMode::per_rank,
                                                     100000,
                                                     10000,
                                                     2,
                                                     "slurm+proc_status"};
    const auto adjusted = ModuleRI::reserve_sternheimer_shared_memory(memory, 10000);
    EXPECT_EQ(adjusted.current_bytes, 20000U);
    EXPECT_EQ(adjusted.limit_bytes, 100000U);
    EXPECT_EQ(ModuleRI::plan_sternheimer_channel_workers(40, 30, 1, 0, adjusted).automatic_workers, 9);
}

TEST(SternheimerChannelResources, ReservesSharedMemoryByReducingAvailableBytes)
{
    const ModuleRI::SternheimerMemorySnapshot memory{ModuleRI::SternheimerMemoryAccountingMode::available,
                                                     100000,
                                                     0,
                                                     2,
                                                     "proc_meminfo"};
    const auto adjusted = ModuleRI::reserve_sternheimer_shared_memory(memory, 10000);
    EXPECT_EQ(adjusted.current_bytes, 0U);
    EXPECT_EQ(adjusted.limit_bytes, 80000U);
    EXPECT_EQ(ModuleRI::plan_sternheimer_channel_workers(40, 30, 1, 0, adjusted).automatic_workers, 15);
    EXPECT_THROW(ModuleRI::reserve_sternheimer_shared_memory(memory, 50001), std::runtime_error);
    const auto exhausted = ModuleRI::reserve_sternheimer_shared_memory(memory, 50000);
    EXPECT_EQ(exhausted.limit_bytes, 0U);
    EXPECT_THROW(ModuleRI::plan_sternheimer_channel_workers(40, 30, 1, 0, exhausted), std::invalid_argument);
}

TEST(SternheimerChannelResources, SharedReservationLeavesUnknownFallbackUnchanged)
{
    ModuleRI::SternheimerMemorySnapshot memory;
    memory.local_mpi_ranks = 3;
    const auto adjusted
        = ModuleRI::reserve_sternheimer_shared_memory(memory, std::numeric_limits<std::uint64_t>::max());
    EXPECT_EQ(adjusted.mode, ModuleRI::SternheimerMemoryAccountingMode::fallback_one);
    EXPECT_EQ(adjusted.limit_bytes, 0U);
    EXPECT_EQ(adjusted.current_bytes, 0U);
    EXPECT_EQ(adjusted.local_mpi_ranks, 3);
    EXPECT_EQ(adjusted.source, "unavailable");
}

TEST(SternheimerChannelResources, ZeroSharedReservationPreservesEveryAccountingMode)
{
    for (const auto mode: {ModuleRI::SternheimerMemoryAccountingMode::node_aggregate,
                          ModuleRI::SternheimerMemoryAccountingMode::per_rank,
                          ModuleRI::SternheimerMemoryAccountingMode::available,
                          ModuleRI::SternheimerMemoryAccountingMode::fallback_one})
    {
        const ModuleRI::SternheimerMemorySnapshot memory{mode, 100000, 10000, 2, "test"};
        const auto adjusted = ModuleRI::reserve_sternheimer_shared_memory(memory, 0);
        EXPECT_EQ(adjusted.mode, memory.mode);
        EXPECT_EQ(adjusted.limit_bytes, memory.limit_bytes);
        EXPECT_EQ(adjusted.current_bytes, memory.current_bytes);
        EXPECT_EQ(adjusted.local_mpi_ranks, memory.local_mpi_ranks);
        EXPECT_EQ(adjusted.source, memory.source);
    }
}

TEST(SternheimerChannelResources, SharedReservationChecksArithmeticOverflowAndRankCount)
{
    const auto maximum = std::numeric_limits<std::uint64_t>::max();
    for (const auto mode: {ModuleRI::SternheimerMemoryAccountingMode::node_aggregate,
                          ModuleRI::SternheimerMemoryAccountingMode::per_rank})
    {
        ModuleRI::SternheimerMemorySnapshot memory{mode, maximum, maximum, 2, "test"};
        EXPECT_THROW(ModuleRI::reserve_sternheimer_shared_memory(memory, 1), std::overflow_error);
        memory.current_bytes = 0;
        if (mode == ModuleRI::SternheimerMemoryAccountingMode::node_aggregate)
        {
            EXPECT_THROW(ModuleRI::reserve_sternheimer_shared_memory(memory, maximum), std::overflow_error);
        }
        else
        {
            EXPECT_EQ(ModuleRI::reserve_sternheimer_shared_memory(memory, maximum).current_bytes, maximum);
        }
        memory.local_mpi_ranks = 0;
        EXPECT_THROW(ModuleRI::reserve_sternheimer_shared_memory(memory, 1), std::invalid_argument);
        memory.local_mpi_ranks = -1;
        EXPECT_THROW(ModuleRI::reserve_sternheimer_shared_memory(memory, 1), std::invalid_argument);
    }
    const ModuleRI::SternheimerMemorySnapshot available{ModuleRI::SternheimerMemoryAccountingMode::available,
                                                       maximum,
                                                       0,
                                                       2,
                                                       "test"};
    EXPECT_THROW(ModuleRI::reserve_sternheimer_shared_memory(available, maximum), std::overflow_error);
}

TEST(SternheimerChannelResources, SharedReservationRejectsPackedProjectorsBeforeAllocation)
{
    constexpr std::uint64_t gib = 1ULL << 30;
    constexpr std::uint64_t shared_bytes = (2ULL * 68 + 676) * 72 * 72 * 135 * sizeof(std::complex<double>);
    const ModuleRI::SternheimerMemorySnapshot memory{ModuleRI::SternheimerMemoryAccountingMode::node_aggregate,
                                                     108 * gib,
                                                     74 * gib,
                                                     1,
                                                     "cgroup_v2"};
    EXPECT_NO_THROW(ModuleRI::plan_sternheimer_channel_workers(3897, 30, 72 * 72 * 135, 1, memory, 2));
    const auto adjusted = ModuleRI::reserve_sternheimer_shared_memory(memory, shared_bytes);
    EXPECT_EQ(adjusted.current_bytes, memory.current_bytes + shared_bytes);
    EXPECT_THROW(ModuleRI::plan_sternheimer_channel_workers(3897, 30, 72 * 72 * 135, 1, adjusted, 2),
                 std::runtime_error);
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
    EXPECT_EQ(parse_sternheimer_kib_field("MemTotal: 4096 kB\nMemAvailable: 2048 kB\n", "MemAvailable"),
              2048ULL * 1024ULL);
    EXPECT_EQ(parse_sternheimer_kib_field("VmPeak:\t2048 kB\nVmRSS:\t1024 kB\n", "VmRSS"), 1024ULL * 1024ULL);
    EXPECT_FALSE(parse_sternheimer_kib_field("VmRSS: invalid kB\n", "VmRSS").has_value());
}

TEST(SternheimerChannelResources, ParsesCgroupMembershipPaths)
{
    EXPECT_EQ(ModuleRI::detail::parse_sternheimer_cgroup_v2_path("0::/slurm/job_7/step_0\n"), "/slurm/job_7/step_0");
    EXPECT_EQ(
        ModuleRI::detail::parse_sternheimer_cgroup_v1_memory_path("4:cpu:/slurm/job_7\n5:memory:/slurm/job_7/step_0\n"),
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
