#include "source_lcao/module_ri/sternheimer_abacus_st_smoke.h"

#include "source_lcao/module_ri/sternheimer_channel_parallel.h"

#include <atomic>
#include <chrono>
#include <complex>
#include <gtest/gtest.h>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

TEST(SternheimerChannelParallel, ExecutesConcurrentlyAndReturnsChannelOrder)
{
#ifdef _OPENMP
    const int previous_threads = omp_get_max_threads();
    const int previous_dynamic = omp_get_dynamic();
    omp_set_dynamic(0);
    omp_set_num_threads(4);
#endif

    std::atomic<int> active_workers{0};
    std::atomic<int> peak_workers{0};
    const std::vector<int> results
        = ModuleRI::run_sternheimer_channel_tasks<int>(32, [&active_workers, &peak_workers](const int channel_index) {
              const int active = active_workers.fetch_add(1) + 1;
              int observed_peak = peak_workers.load();
              while (active > observed_peak && !peak_workers.compare_exchange_weak(observed_peak, active))
              {
              }
              std::this_thread::sleep_for(std::chrono::milliseconds(2));
              active_workers.fetch_sub(1);
              return channel_index * channel_index;
          });

    ASSERT_EQ(results.size(), 32U);
    for (int channel_index = 0; channel_index != 32; ++channel_index)
    {
        EXPECT_EQ(results[static_cast<std::size_t>(channel_index)], channel_index * channel_index);
    }
#ifdef _OPENMP
    EXPECT_GT(peak_workers.load(), 1);
    omp_set_num_threads(previous_threads);
    omp_set_dynamic(previous_dynamic);
#else
    EXPECT_EQ(peak_workers.load(), 1);
#endif
}

TEST(SternheimerChannelParallel, WritesIndependentColumnsLikeSerialExecution)
{
    constexpr int num_channels = 17;
    std::vector<double> serial(num_channels * num_channels, 0.0);
    std::vector<double> parallel(num_channels * num_channels, 0.0);
    for (int column = 0; column != num_channels; ++column)
    {
        for (int row = 0; row != num_channels; ++row)
        {
            serial[static_cast<std::size_t>(row) * num_channels + column] = 1000.0 * row + column;
        }
    }

    static_cast<void>(ModuleRI::run_sternheimer_channel_tasks<int>(
        num_channels,
        [&parallel](const int column) {
            for (int row = 0; row != num_channels; ++row)
            {
                parallel[static_cast<std::size_t>(row) * num_channels + column] = 1000.0 * row + column;
            }
            return column;
        }));
    EXPECT_EQ(parallel, serial);
}

TEST(SternheimerChannelParallel, RethrowsFirstIndexedExceptionAfterAllTasksFinish)
{
    std::atomic<int> completed_tasks{0};
    try
    {
        static_cast<void>(ModuleRI::run_sternheimer_channel_tasks<int>(16, [&completed_tasks](const int channel_index) {
            completed_tasks.fetch_add(1);
            if (channel_index == 3 || channel_index == 9)
            {
                throw std::runtime_error("channel " + std::to_string(channel_index));
            }
            return channel_index;
        }));
        FAIL() << "Expected a channel task exception.";
    }
    catch (const std::runtime_error& error)
    {
        EXPECT_STREQ(error.what(), "channel 3");
    }
    EXPECT_EQ(completed_tasks.load(), 16);
    EXPECT_THROW(ModuleRI::run_sternheimer_channel_tasks<int>(-1, [](const int) { return 0; }), std::invalid_argument);
}

TEST(SternheimerABACUSSTSmoke, FormatsLinearResponseReport)
{
    ModuleRI::SternheimerABACUSSTSmokeResult result;
    result.grid_data.grid.nx = 2;
    result.grid_data.grid.ny = 3;
    result.grid_data.grid.nz = 4;
    result.grid_data.volume_element = 0.125;
    result.omega = 0.5;
    result.pca_threshold = 1.0e-4;
    result.ccp_rmesh_times = 10.0;
    result.perturbation_source = "abfs_ccp";
    result.num_available_channels = 5;

    ModuleRI::SternheimerABACUSSTChannelResult channel;
    channel.band_index = 0;
    channel.channel_index = 1;
    channel.atom_index = 0;
    channel.angular_momentum = 1;
    channel.radial_index = 2;
    channel.magnetic_index = 1;
    channel.fd_eigenvalue = -1.25;
    channel.occupation = 2.0;
    channel.rhs_norm = 0.4;
    channel.projected_rhs_norm = 0.3;
    channel.solver_converged = true;
    channel.solver_iterations = 7;
    channel.solver_relative_residual = 1.0e-9;
    channel.equation_residual_norm = 2.0e-9;
    channel.polarizability = std::complex<double>(-0.12, 0.03);
    result.channels.push_back(channel);

    const std::string report = ModuleRI::format_sternheimer_abacus_st_report(result);

    EXPECT_NE(report.find("# ABACUS Sternheimer FD linear-response smoke test"), std::string::npos);
    EXPECT_NE(report.find("grid 2 3 4 size 24 dV 0.125"), std::string::npos);
    EXPECT_NE(report.find("omega_Ry 0.5"), std::string::npos);
    EXPECT_NE(report.find("pca_threshold 0.0001"), std::string::npos);
    EXPECT_NE(report.find("perturbation_source abfs_ccp"), std::string::npos);
    EXPECT_NE(report.find("available_channels 5"), std::string::npos);
    EXPECT_NE(report.find("0 1 0 1 2 1 -1.25 2 0.4 0.3 yes 7 1e-09 2e-09 -0.12 0.03"),
              std::string::npos);
}

namespace
{

ModuleRI::SternheimerLCAOOccupiedKPoint make_occupied_kpoint(const int local_k_index,
                                                             const int global_k_index,
                                                             const int spin_index,
                                                             const ModuleRI::SternheimerReducedKPoint& kpoint,
                                                             const double kweight)
{
    ModuleRI::SternheimerLCAOOccupiedKPoint record;
    record.local_k_index = local_k_index;
    record.global_k_index = global_k_index;
    record.spin_index = spin_index;
    record.kpoint = kpoint;
    record.kweight = kweight;
    record.eigenvalues = {-1.0};
    record.occupations = {1.0};
    record.coefficients = {{std::complex<double>(1.0, 0.0),
                            std::complex<double>(0.0, 0.0),
                            std::complex<double>(0.0, 0.0)}};
    return record;
}

} // namespace

TEST(SternheimerABACUSSTSmoke, LimitsDiagnosticChannelsPerAtomWithoutEmptyBlocks)
{
    std::vector<ModuleRI::SternheimerABFBlochGridChannel> channels(4);
    channels[0].atom_index = 0;
    channels[0].atom_local_index = 0;
    channels[1].atom_index = 0;
    channels[1].atom_local_index = 1;
    channels[2].atom_index = 1;
    channels[2].atom_local_index = 0;
    channels[3].atom_index = 1;
    channels[3].atom_local_index = 1;

    const auto limited = ModuleRI::limit_sternheimer_abf_channels_per_atom(channels, 1);

    ASSERT_EQ(limited.size(), 2);
    EXPECT_EQ(limited[0].atom_index, 0);
    EXPECT_EQ(limited[0].atom_local_index, 0);
    EXPECT_EQ(limited[0].channel_index, 0);
    EXPECT_EQ(limited[1].atom_index, 1);
    EXPECT_EQ(limited[1].atom_local_index, 0);
    EXPECT_EQ(limited[1].channel_index, 1);
}

TEST(SternheimerABACUSSTSmoke, DistinguishesTwoKPointsFromTwoSpinChannels)
{
    const auto k0 = make_occupied_kpoint(0, 0, 0, {0.0, 0.0, 0.0}, 0.5);
    const auto k1 = make_occupied_kpoint(1, 1, 0, {0.5, 0.0, 0.0}, 0.5);
    const std::vector<ModuleRI::SternheimerLCAOOccupiedKPoint> two_kpoints = {k0, k1};
    EXPECT_NO_THROW(ModuleRI::validate_sternheimer_lcao_occupied_kpoints(two_kpoints, 2, 2, 1, 3));
    EXPECT_EQ(ModuleRI::sternheimer_lcao_total_occupied_bands(two_kpoints), 2);

    auto spin_down = k0;
    spin_down.local_k_index = 1;
    spin_down.global_k_index = 1;
    spin_down.spin_index = 1;
    const std::vector<ModuleRI::SternheimerLCAOOccupiedKPoint> two_spins = {k0, spin_down};
    EXPECT_NO_THROW(ModuleRI::validate_sternheimer_lcao_occupied_kpoints(two_spins, 2, 2, 2, 3));
}

TEST(SternheimerABACUSSTSmoke, SelectsEveryGammaSpinRecordForMolecularResponse)
{
    const auto spin_up = make_occupied_kpoint(0, 0, 0, {0.0, 0.0, 0.0}, 1.0);
    auto spin_down = make_occupied_kpoint(1, 1, 1, {0.0, 0.0, 0.0}, 1.0);

    const std::vector<ModuleRI::SternheimerLCAOOccupiedKPoint> records = {spin_down, spin_up};
    const auto selected = ModuleRI::select_sternheimer_gamma_spin_records(records, 2);

    ASSERT_EQ(selected.size(), 2U);
    EXPECT_EQ(selected[0]->spin_index, 0);
    EXPECT_EQ(selected[0]->local_k_index, 0);
    EXPECT_EQ(selected[1]->spin_index, 1);
    EXPECT_EQ(selected[1]->local_k_index, 1);

    spin_down.spin_index = 0;
    EXPECT_THROW(ModuleRI::select_sternheimer_gamma_spin_records({spin_up, spin_down}, 2),
                 std::invalid_argument);
}

TEST(SternheimerABACUSSTSmoke, RejectsDuplicateOrIncompleteGlobalKRecords)
{
    const auto k0 = make_occupied_kpoint(0, 0, 0, {0.0, 0.0, 0.0}, 0.5);
    const auto k1 = make_occupied_kpoint(1, 1, 0, {0.5, 0.0, 0.0}, 0.5);

    auto duplicate = k1;
    duplicate.global_k_index = 0;
    EXPECT_THROW(ModuleRI::validate_sternheimer_lcao_occupied_kpoints({k0, duplicate}, 2, 2, 1, 3),
                 std::invalid_argument);
    EXPECT_THROW(ModuleRI::validate_sternheimer_lcao_occupied_kpoints({k0}, 2, 2, 1, 3),
                 std::invalid_argument);
}

TEST(SternheimerABACUSSTSmoke, GammaRecordPreservesLegacyWeightedOccupation)
{
    const auto gamma = make_occupied_kpoint(0, 0, 0, {0.0, 0.0, 0.0}, 2.0);
    EXPECT_NO_THROW(ModuleRI::validate_sternheimer_lcao_occupied_kpoints({gamma}, 1, 1, 1, 3));
    EXPECT_DOUBLE_EQ(ModuleRI::sternheimer_lcao_weighted_occupation(gamma, 0), 2.0);
}

TEST(SternheimerABACUSSTSmoke, ValidatesOptionalUnoccupiedLCAOStates)
{
    auto record = make_occupied_kpoint(0, 0, 0, {0.0, 0.0, 0.0}, 2.0);
    record.unoccupied_eigenvalues = {0.5};
    record.unoccupied_coefficients = {{std::complex<double>(0.0, 0.0),
                                       std::complex<double>(1.0, 0.0),
                                       std::complex<double>(0.0, 0.0)}};
    EXPECT_NO_THROW(ModuleRI::validate_sternheimer_lcao_occupied_kpoints({record}, 1, 1, 1, 3));

    record.unoccupied_coefficients[0].pop_back();
    EXPECT_THROW(ModuleRI::validate_sternheimer_lcao_occupied_kpoints({record}, 1, 1, 1, 3),
                 std::invalid_argument);
}

TEST(SternheimerABACUSSTSmoke, BuildsTwoKPointNonzeroQResponsePlan)
{
    const auto k0 = make_occupied_kpoint(0, 0, 0, {0.0, 0.0, 0.0}, 1.0);
    const auto k1 = make_occupied_kpoint(1, 1, 0, {0.5, 0.0, 0.0}, 1.0);

    const auto plan = ModuleRI::build_sternheimer_periodic_response_plan({k1, k0}, 2);

    EXPECT_EQ(plan.iq, 2);
    EXPECT_EQ(plan.qpoint, (ModuleRI::SternheimerReducedKPoint{0.5, 0.0, 0.0}));
    EXPECT_EQ(plan.record_index_by_global_k, (std::vector<int>{1, 0}));
    ASSERT_EQ(plan.kq_pairs.size(), 2);
    EXPECT_EQ(plan.kq_pairs[0].source_index, 0);
    EXPECT_EQ(plan.kq_pairs[0].target_index, 1);
    EXPECT_EQ(plan.kq_pairs[1].source_index, 1);
    EXPECT_EQ(plan.kq_pairs[1].target_index, 0);
    EXPECT_DOUBLE_EQ(plan.kweight_sum, 2.0);
}

TEST(SternheimerABACUSSTSmoke, BuildsPeriodicGammaResponsePlanBySelfMappingKPoints)
{
    auto k0 = make_occupied_kpoint(0, 0, 0, {0.0, 0.0, 0.0}, 1.0);
    const auto k1 = make_occupied_kpoint(1, 1, 0, {0.5, 0.0, 0.0}, 1.0);

    const auto plan = ModuleRI::build_sternheimer_periodic_response_plan({k0, k1}, 1);

    EXPECT_EQ(plan.iq, 1);
    EXPECT_EQ(plan.qpoint, (ModuleRI::SternheimerReducedKPoint{0.0, 0.0, 0.0}));
    ASSERT_EQ(plan.kq_pairs.size(), 2);
    EXPECT_EQ(plan.kq_pairs[0].source_index, 0);
    EXPECT_EQ(plan.kq_pairs[0].target_index, 0);
    EXPECT_EQ(plan.kq_pairs[1].source_index, 1);
    EXPECT_EQ(plan.kq_pairs[1].target_index, 1);
}

TEST(SternheimerABACUSSTSmoke, RejectsFractionalOccupationForPeriodicResponse)
{
    auto k0 = make_occupied_kpoint(0, 0, 0, {0.0, 0.0, 0.0}, 1.0);
    const auto k1 = make_occupied_kpoint(1, 1, 0, {0.5, 0.0, 0.0}, 1.0);

    k0.occupations[0] = 0.5;
    EXPECT_THROW(ModuleRI::build_sternheimer_periodic_response_plan({k0, k1}, 2), std::invalid_argument);
}

TEST(SternheimerABACUSSTSmoke, ZeroQIndexPreservesSingleGammaPlan)
{
    const auto gamma = make_occupied_kpoint(0, 0, 0, {0.0, 0.0, 0.0}, 2.0);
    const auto plan = ModuleRI::build_sternheimer_periodic_response_plan({gamma}, 0);

    EXPECT_EQ(plan.iq, 1);
    EXPECT_EQ(plan.qpoint, (ModuleRI::SternheimerReducedKPoint{0.0, 0.0, 0.0}));
    ASSERT_EQ(plan.kq_pairs.size(), 1);
    EXPECT_EQ(plan.kq_pairs[0].source_index, 0);
    EXPECT_EQ(plan.kq_pairs[0].target_index, 0);
    EXPECT_DOUBLE_EQ(plan.kweight_sum, 2.0);
}

TEST(SternheimerABACUSSTSmoke, SelectsMassiddaFactorOnlyForPeriodicGamma)
{
    constexpr double massidda_chi = 1.25;
    EXPECT_DOUBLE_EQ(ModuleRI::sternheimer_periodic_gamma_inverse_k2(
                         {0.0, 0.0, 0.0}, "massidda", massidda_chi),
                     massidda_chi);
    EXPECT_DOUBLE_EQ(ModuleRI::sternheimer_periodic_gamma_inverse_k2(
                         {0.25, 0.0, 0.0}, "massidda", massidda_chi),
                     0.0);
    EXPECT_THROW(ModuleRI::sternheimer_periodic_gamma_inverse_k2(
                     {0.0, 0.0, 0.0}, "limits", massidda_chi),
                 std::invalid_argument);
    EXPECT_THROW(ModuleRI::sternheimer_periodic_gamma_inverse_k2(
                     {0.0, 0.0, 0.0}, "massidda", 0.0),
                 std::invalid_argument);
}

TEST(SternheimerABACUSSTSmoke, ValidatesPeriodicMonkhorstPackDimensions)
{
    EXPECT_NO_THROW(ModuleRI::validate_sternheimer_periodic_kmesh({4, 4, 4}, 64));
    EXPECT_THROW(ModuleRI::validate_sternheimer_periodic_kmesh({4, 4, 4}, 63),
                 std::invalid_argument);
    EXPECT_THROW(ModuleRI::validate_sternheimer_periodic_kmesh({4, 0, 4}, 64),
                 std::invalid_argument);
}

TEST(SternheimerABACUSSTSmoke, AssignsContiguousKPointOwners)
{
    const int global_kpoint_count = 10;
    const int kpoint_groups = 4;
    const std::vector<int> expected_owners = {0, 0, 0, 1, 1, 1, 2, 2, 3, 3};

    for (int ik = 0; ik != global_kpoint_count; ++ik)
    {
        EXPECT_EQ(ModuleRI::sternheimer_kpoint_owner_group(ik, global_kpoint_count, kpoint_groups),
                  expected_owners[static_cast<std::size_t>(ik)]);
    }
    EXPECT_THROW(ModuleRI::sternheimer_kpoint_owner_group(-1, global_kpoint_count, kpoint_groups),
                 std::invalid_argument);
    EXPECT_THROW(ModuleRI::sternheimer_kpoint_owner_group(0, global_kpoint_count, 0),
                 std::invalid_argument);
    EXPECT_THROW(ModuleRI::sternheimer_kpoint_owner_group(0, 3, 4), std::invalid_argument);
}

TEST(SternheimerABACUSSTSmoke, PartitionsResponsePairsBySourceKWithoutOverlap)
{
    std::vector<ModuleRI::SternheimerLCAOOccupiedKPoint> records;
    for (int ik = 0; ik != 8; ++ik)
    {
        records.push_back(make_occupied_kpoint(ik,
                                               ik,
                                               0,
                                               {static_cast<double>(ik) / 8.0, 0.0, 0.0},
                                               0.25));
    }
    const auto plan = ModuleRI::build_sternheimer_periodic_response_plan(records, 2);

    std::vector<int> seen(plan.kq_pairs.size(), 0);
    for (int group = 0; group != 4; ++group)
    {
        const auto owned = ModuleRI::sternheimer_owned_kq_pair_indices(plan, group, 4);
        ASSERT_EQ(owned.size(), 2);
        for (const std::size_t pair_index: owned)
        {
            ASSERT_LT(pair_index, plan.kq_pairs.size());
            ++seen[pair_index];
            EXPECT_EQ(ModuleRI::sternheimer_kpoint_owner_group(plan.kq_pairs[pair_index].source_index, 8, 4),
                      group);
        }
    }
    EXPECT_EQ(seen, std::vector<int>(plan.kq_pairs.size(), 1));
}
