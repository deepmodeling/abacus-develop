#include "source_lcao/module_ri/sternheimer_abacus_st_smoke.h"

#include <gtest/gtest.h>
#include <algorithm>
#include <numeric>
#include <set>

namespace
{
using ScheduleUnderTest = ModuleRI::SternheimerKPointSchedule;

ModuleRI::SternheimerPeriodicResponsePlan make_plan(const int iq)
{
    std::vector<ModuleRI::SternheimerLCAOOccupiedKPoint> records(16);
    for (int ik = 0; ik != 16; ++ik)
    {
        auto& record = records[ik];
        record.local_k_index = record.global_k_index = record.zero_order_k_index = ik;
        record.spin_index = 0;
        record.kpoint = {static_cast<double>(ik / 4) / 4.0, static_cast<double>(ik % 4) / 4.0, 0.0};
        record.kweight = 1.0 / 16.0;
        record.occupations = {1.0};
    }
    return ModuleRI::build_sternheimer_periodic_response_plan(records, iq);
}

std::vector<bool> tr_representatives()
{
    std::vector<int> identity(16), reversal(16);
    std::iota(identity.begin(), identity.end(), 0);
    for (int ik = 0; ik != 16; ++ik)
    {
        reversal[ik] = ((4 - ik / 4) % 4) * 4 + (4 - ik % 4) % 4;
    }
    const auto orbits = ModuleRI::build_sternheimer_fixed_q_k_orbits_from_permutations(
        16, {identity, reversal});
    std::vector<bool> active(16, false);
    for (const auto& orbit: orbits)
    {
        active[orbit.representative_ik_full] = true;
    }
    return active;
}
} // namespace

TEST(SternheimerRepresentativeOwnership, BalancesTenRepresentativesAcrossFiveGroups)
{
    RecordProperty("tr_q_indices", "1,3,9,11");
    RecordProperty("representative_count", 10);
    RecordProperty("group_count", 5);
    const auto active = tr_representatives();
    std::vector<int> representatives;
    for (int ik = 0; ik != 16; ++ik)
    {
        if (active[ik]) representatives.push_back(ik + 1);
    }
    EXPECT_EQ(representatives, (std::vector<int>{1, 2, 3, 5, 6, 7, 8, 9, 10, 11}));
    const ScheduleUnderTest schedule(active, 5);
    for (const int iq: {1, 3, 9, 11})
    {
        auto plan = make_plan(iq);
        std::reverse(plan.kq_pairs.begin(), plan.kq_pairs.end());
        const auto original = plan;
        std::vector<int> seen(16, 0);
        for (int group = 0; group != 5; ++group)
        {
            SCOPED_TRACE(::testing::Message() << "q=" << iq << " group=" << group);
            const auto owned = schedule.owned_pair_indices(plan, group);
            EXPECT_EQ(owned.size(), 2U);
            for (const auto pair_index: owned)
            {
                ASSERT_LT(pair_index, plan.kq_pairs.size());
                const auto& pair = plan.kq_pairs[pair_index];
                EXPECT_TRUE(active[pair.source_index]);
                EXPECT_EQ(schedule.owner_group(pair.source_index), group);
                EXPECT_EQ(pair.source_index, original.kq_pairs[pair_index].source_index);
                EXPECT_EQ(pair.target_index, original.kq_pairs[pair_index].target_index);
                EXPECT_EQ(pair.reciprocal_shift, original.kq_pairs[pair_index].reciprocal_shift);
                ++seen[pair.source_index];
            }
        }
        for (int ik = 0; ik != 16; ++ik) EXPECT_EQ(seen[ik], active[ik] ? 1 : 0);
        EXPECT_EQ(plan.iq, iq);
        EXPECT_EQ(plan.kweight_sum, original.kweight_sum);
    }
}

TEST(SternheimerRepresentativeOwnership, OwnsEveryPairFrequencyReplicaOnceWithTwoPairsPerRank)
{
    RecordProperty("base_nproc", 30);
    RecordProperty("nfreq", 6);
    RecordProperty("pairs_per_rank", 2);
    RecordProperty("unique_fragments_per_q", 60);
    RecordProperty("replica_counts_tested", "1,2,3");
    const auto active = tr_representatives();
    const ScheduleUnderTest schedule(active, 5);
    constexpr int nfreq = 6;
    for (const int iq: {1, 3, 9, 11})
    {
        const auto plan = make_plan(iq);
        for (const int replicas: {1, 2, 3})
        {
            for (const int shift: {-7, 0, 1, 7})
            {
                const int ranks = 5 * nfreq * replicas;
                std::vector<int> seen(16 * nfreq * replicas, 0), writers(16 * nfreq, 0);
                for (int rank = 0; rank != ranks; ++rank)
                {
                    const auto layout = ModuleRI::sternheimer_nested_mpi_replica_layout(
                        5, nfreq, ranks, rank, replicas > 1);
                    const int group = layout.local_response_slot / nfreq;
                    const auto owned = schedule.owned_pair_indices(plan, group);
                    int rank_pairs = 0;
                    for (const auto pair_index: owned)
                    {
                        const int source = plan.kq_pairs[pair_index].source_index;
                        for (int frequency = 0; frequency != nfreq; ++frequency)
                        {
                            const auto assignment = schedule.assignment(source, frequency, nfreq, 5 * nfreq, shift);
                            const int leader = assignment.owner_rank * replicas;
                            if (rank < leader || rank >= leader + replicas) continue;
                            EXPECT_EQ(assignment.kpoint_group, group);
                            EXPECT_EQ(assignment.frequency_slot, layout.local_response_slot % nfreq);
                            ++seen[(source * nfreq + frequency) * replicas + layout.local_replica];
                            if (layout.local_replica == 0) ++writers[source * nfreq + frequency];
                            ++rank_pairs;
                        }
                    }
                    EXPECT_EQ(rank_pairs, 2) << "q=" << iq << " rank=" << rank
                                             << " replicas=" << replicas << " shift=" << shift;
                }
                for (int source = 0; source != 16; ++source)
                {
                    for (int frequency = 0; frequency != nfreq; ++frequency)
                    {
                        EXPECT_EQ(writers[source * nfreq + frequency], active[source] ? 1 : 0);
                        for (int replica = 0; replica != replicas; ++replica)
                            EXPECT_EQ(seen[(source * nfreq + frequency) * replicas + replica], active[source] ? 1 : 0);
                    }
                }
            }
        }
    }
}

TEST(SternheimerRepresentativeOwnership, PreservesFullGridPairAndFrequencyOwners)
{
    RecordProperty("full_k_count", 16);
    RecordProperty("group_count", 8);
    RecordProperty("nontr_mapping", "unchanged");
    const ScheduleUnderTest schedule(std::vector<bool>(16, true), 8);
    for (const int iq: {5, 6, 7, 8, 10})
    {
        const auto plan = make_plan(iq);
        for (int group = 0; group != 8; ++group)
        {
            const auto owned = schedule.owned_pair_indices(plan, group);
            EXPECT_EQ(owned, ModuleRI::sternheimer_owned_kq_pair_indices(plan, group, 8));
            ASSERT_EQ(owned.size(), 2U);
            for (const auto pair_index: owned)
            {
                const int source = plan.kq_pairs[pair_index].source_index;
                EXPECT_EQ(schedule.owner_group(source), ModuleRI::sternheimer_kpoint_owner_group(source, 16, 8));
                for (int frequency = 0; frequency != 6; ++frequency)
                {
                    for (const int shift: {-7, 0, 1, 7})
                    {
                        const auto actual = schedule.assignment(source, frequency, 6, 48, shift);
                        const auto original = ModuleRI::sternheimer_nested_mpi_assignment(source, 16, frequency, 6, 8, 48, shift);
                        EXPECT_EQ(actual.kpoint_group, original.kpoint_group);
                        EXPECT_EQ(actual.frequency_slot, original.frequency_slot);
                        EXPECT_EQ(actual.owner_rank, original.owner_rank);
                    }
                }
            }
        }
    }
}

TEST(SternheimerRepresentativeOwnership, RejectsIdleGroupsAndInactiveOwners)
{
    const auto active = tr_representatives();
    EXPECT_THROW(ScheduleUnderTest(active, 0), std::invalid_argument);
    EXPECT_THROW(ScheduleUnderTest(active, -1), std::invalid_argument);
    EXPECT_THROW(ScheduleUnderTest(active, 11), std::invalid_argument);
    EXPECT_THROW(ScheduleUnderTest(std::vector<bool>(16, false), 1), std::invalid_argument);
    EXPECT_THROW(ScheduleUnderTest(std::vector<bool>(), 1), std::invalid_argument);
    const ScheduleUnderTest schedule(active, 5);
    EXPECT_THROW(schedule.owner_group(-1), std::invalid_argument);
    EXPECT_THROW(schedule.owner_group(16), std::invalid_argument);
    EXPECT_THROW(schedule.owner_group(3), std::invalid_argument);
    EXPECT_THROW(schedule.assignment(3, 0, 6, 30), std::invalid_argument);
    EXPECT_THROW(schedule.assignment(0, 0, 6, 48), std::invalid_argument);
    EXPECT_THROW(schedule.assignment(0, 6, 6, 30), std::invalid_argument);
    const auto plan = make_plan(3);
    EXPECT_THROW(schedule.owned_pair_indices(plan, -1), std::invalid_argument);
    EXPECT_THROW(schedule.owned_pair_indices(plan, 5), std::invalid_argument);
    auto malformed = plan;
    malformed.record_index_by_global_k.pop_back();
    EXPECT_THROW(schedule.owned_pair_indices(malformed, 0), std::invalid_argument);
    malformed = plan;
    malformed.kq_pairs.front().source_index = 16;
    EXPECT_THROW(schedule.owned_pair_indices(malformed, 0), std::invalid_argument);
}

TEST(SternheimerRepresentativeOwnership, CoversChannelBatchesExactlyOnceAcrossReplicas)
{
    const auto active = tr_representatives();
    const ScheduleUnderTest schedule(active, 5);
    const auto plan = make_plan(3);
    constexpr int nfreq = 6;
    constexpr int occupied_states = 3;
    constexpr int batches = 7;
    for (const int replicas: {1, 2, 4})
    {
        const int ranks = 30 * replicas;
        std::vector<int> seen(16 * nfreq * occupied_states * batches, 0);
        for (int rank = 0; rank != ranks; ++rank)
        {
            const auto layout = ModuleRI::sternheimer_nested_mpi_replica_layout(5, nfreq, ranks, rank, true);
            const auto owned = schedule.owned_pair_indices(plan, layout.local_response_slot / nfreq);
            for (const auto pair_index: owned)
            {
                const int source = plan.kq_pairs[pair_index].source_index;
                for (int frequency = 0; frequency != nfreq; ++frequency)
                {
                    const auto assignment = schedule.assignment(source, frequency, nfreq, 30, 1);
                    if (layout.local_response_slot != assignment.owner_rank) continue;
                    for (int occupied = 0; occupied != occupied_states; ++occupied)
                    {
                        for (int batch = 0; batch != batches; ++batch)
                        {
                            if (ModuleRI::sternheimer_channel_batch_replica_owner(occupied, batch, batches, replicas)
                                != layout.local_replica) continue;
                            ++seen[((source * nfreq + frequency) * occupied_states + occupied) * batches + batch];
                        }
                    }
                }
            }
        }
        for (int source = 0; source != 16; ++source)
        {
            for (int task = 0; task != nfreq * occupied_states * batches; ++task)
                EXPECT_EQ(seen[source * nfreq * occupied_states * batches + task], active[source] ? 1 : 0);
        }
    }
}

TEST(SternheimerRepresentativeOwnership, HandlesSparseNonzeroFirstRepresentativeAndSingleGroup)
{
    std::vector<bool> active(16, false);
    for (const int source: {2, 7, 12}) active[source] = true;
    const auto plan = make_plan(9);
    const ScheduleUnderTest schedule(active, 3);
    for (int group = 0; group != 3; ++group)
    {
        const auto owned = schedule.owned_pair_indices(plan, group);
        ASSERT_EQ(owned.size(), 1U);
        EXPECT_EQ(plan.kq_pairs[owned.front()].source_index, (std::vector<int>{2, 7, 12})[group]);
        EXPECT_EQ(schedule.assignment(plan.kq_pairs[owned.front()].source_index, 0, 6, 18).owner_rank, group * 6);
    }
    const ScheduleUnderTest single_group(active, 1);
    EXPECT_EQ(single_group.owned_pair_indices(plan, 0).size(), 3U);
    for (const int source: {2, 7, 12}) EXPECT_EQ(single_group.owner_group(source), 0);
}

TEST(SternheimerRepresentativeOwnership, WritesSixtyUniqueFragmentsWithOriginalGlobalIndices)
{
    RecordProperty("output_index_space", "original_global_k");
    const auto active = tr_representatives();
    const ScheduleUnderTest schedule(active, 5);
    for (const int iq: {1, 3, 9, 11})
    {
        const auto plan = make_plan(iq);
        std::set<std::string> filenames;
        for (int rank = 0; rank != 30; ++rank)
        {
            const auto owned = schedule.owned_pair_indices(plan, rank / 6);
            for (const auto pair_index: owned)
            {
                const int source = plan.kq_pairs[pair_index].source_index;
                for (int frequency = 0; frequency != 6; ++frequency)
                {
                    if (schedule.assignment(source, frequency, 6, 30, 1).owner_rank != rank) continue;
                    const auto record = ModuleRI::make_sternheimer_partial_response_record(
                        iq, source, frequency + 1, {std::complex<double>(1.0, 2.0)}, 1);
                    EXPECT_EQ(record.iq, iq);
                    EXPECT_EQ(record.ik_full, source);
                    EXPECT_EQ(record.ifrequency, frequency + 1);
                    EXPECT_EQ(record.matrix, (std::vector<std::complex<double>>{{2.0, 0.0}}));
                    EXPECT_TRUE(filenames.insert(record.filename).second);
                }
            }
        }
        EXPECT_EQ(filenames.size(), 60U);
        for (int source = 0; source != 16; ++source)
        {
            for (int frequency = 1; frequency != 7; ++frequency)
                EXPECT_EQ(filenames.count(ModuleRI::sternheimer_partial_response_filename(iq, source, frequency)),
                          active[source] ? 1U : 0U);
        }
    }
}
