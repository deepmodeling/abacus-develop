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

TEST(SternheimerChannelParallel, HonorsExplicitSingleWorkerLimit)
{
#ifdef _OPENMP
    const int previous_threads = omp_get_max_threads();
    const int previous_dynamic = omp_get_dynamic();
    omp_set_dynamic(0);
    omp_set_num_threads(4);
#endif

    std::atomic<int> active_workers{0};
    std::atomic<int> peak_workers{0};
    static_cast<void>(ModuleRI::run_sternheimer_channel_tasks<int>(
        16,
        [&active_workers, &peak_workers](const int channel_index) {
            const int active = active_workers.fetch_add(1) + 1;
            int observed_peak = peak_workers.load();
            while (active > observed_peak && !peak_workers.compare_exchange_weak(observed_peak, active))
            {
            }
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
            active_workers.fetch_sub(1);
            return channel_index;
        },
        1));
    EXPECT_EQ(peak_workers.load(), 1);
#ifdef _OPENMP
    omp_set_num_threads(previous_threads);
    omp_set_dynamic(previous_dynamic);
#endif
}

TEST(SternheimerChannelParallel, SingleOuterWorkerLeavesGridParallelismAvailable)
{
#ifdef _OPENMP
    const int previous_threads = omp_get_max_threads();
    const int previous_dynamic = omp_get_dynamic();
    const int previous_active_levels = omp_get_max_active_levels();
    omp_set_dynamic(0);
    omp_set_num_threads(4);
    omp_set_max_active_levels(1);

    int inner_threads = 0;
    static_cast<void>(ModuleRI::run_sternheimer_channel_tasks<int>(
        1,
        [&inner_threads](const int channel_index) {
#pragma omp parallel num_threads(4)
            {
#pragma omp single
                inner_threads = omp_get_num_threads();
            }
            return channel_index;
        },
        1));
    EXPECT_EQ(inner_threads, 4);

    omp_set_max_active_levels(previous_active_levels);
    omp_set_num_threads(previous_threads);
    omp_set_dynamic(previous_dynamic);
#endif
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
    EXPECT_THROW(ModuleRI::run_sternheimer_channel_tasks<int>(1, [](const int) { return 0; }, -1),
                 std::invalid_argument);
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
    record.zero_order_k_index = local_k_index;
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
    spin_down.zero_order_k_index = 1;
    spin_down.spin_index = 1;
    const std::vector<ModuleRI::SternheimerLCAOOccupiedKPoint> two_spins = {k0, spin_down};
    EXPECT_NO_THROW(ModuleRI::validate_sternheimer_lcao_occupied_kpoints(two_spins, 2, 2, 2, 3));
}

TEST(SternheimerABACUSSTSmoke, AllowsFullKRecordsToShareAnIBZZeroOrderState)
{
    auto k0 = make_occupied_kpoint(0, 0, 0, {0.0, 0.0, 0.0}, 0.5);
    auto k1 = make_occupied_kpoint(1, 1, 0, {0.5, 0.0, 0.0}, 0.5);
    k0.zero_order_k_index = 0;
    k1.zero_order_k_index = 0;

    EXPECT_NO_THROW(
        ModuleRI::validate_sternheimer_lcao_occupied_kpoints({k0, k1}, 2, 2, 1, 3, 1));
    k1.zero_order_k_index = 1;
    EXPECT_THROW(ModuleRI::validate_sternheimer_lcao_occupied_kpoints({k0, k1}, 2, 2, 1, 3, 1),
                 std::invalid_argument);
}

TEST(SternheimerABACUSSTSmoke, ValidatesFullKMeshAgainstSeparateIBZZeroOrderMesh)
{
    constexpr int full_kpoint_count = 8;
    constexpr int ibz_kpoint_count = 3;
    const std::array<int, full_kpoint_count> ibz_index_by_full_k = {0, 1, 1, 2, 1, 2, 2, 1};
    std::vector<ModuleRI::SternheimerLCAOOccupiedKPoint> full_records;
    full_records.reserve(full_kpoint_count);
    for (int full_k_index = 0; full_k_index != full_kpoint_count; ++full_k_index)
    {
        auto record = make_occupied_kpoint(full_k_index,
                                           full_k_index,
                                           0,
                                           {0.0, 0.0, 0.0},
                                           2.0 / static_cast<double>(full_kpoint_count));
        record.zero_order_k_index = ibz_index_by_full_k[static_cast<std::size_t>(full_k_index)];
        full_records.push_back(std::move(record));
    }

    EXPECT_NO_THROW(ModuleRI::validate_sternheimer_full_lcao_occupied_kpoints(
        full_records, ibz_kpoint_count, 1, 3));
}

TEST(SternheimerABACUSSTSmoke, FullKRecordKeepsNormalizedCoefficientsAndUsesUniformWeight)
{
    auto ibz = make_occupied_kpoint(0, 0, 0, {0.0, 0.0, 0.0}, 1.5);
    ibz.zero_order_k_index = 0;
    const std::vector<std::vector<std::complex<double>>> rotated = {
        {{0.0, 1.0}, {1.0, 0.0}, {0.0, 0.0}}};

    const auto full = ModuleRI::make_sternheimer_full_kpoint_record(
        ibz, 3, {0.25, 0.25, 0.0}, 0.25, rotated);

    EXPECT_EQ(full.local_k_index, 3);
    EXPECT_EQ(full.global_k_index, 3);
    EXPECT_EQ(full.zero_order_k_index, 0);
    EXPECT_EQ(full.kpoint, (ModuleRI::SternheimerReducedKPoint{0.25, 0.25, 0.0}));
    EXPECT_DOUBLE_EQ(full.kweight, 0.25);
    EXPECT_EQ(full.coefficients, rotated);
}

TEST(SternheimerABACUSSTSmoke, ListsOneCanonicalFullQPointPerZeroOrderStar)
{
    auto q0 = make_occupied_kpoint(0, 0, 0, {0.0, 0.0, 0.0}, 1.0);
    auto q1 = make_occupied_kpoint(1, 1, 0, {0.5, 0.0, 0.0}, 1.0);
    auto q1_partner = make_occupied_kpoint(2, 2, 0, {-0.5, 0.0, 0.0}, 1.0);
    q1_partner.zero_order_k_index = 1;
    q1_partner.symmetry_spatial_isym = 3;

    EXPECT_EQ(ModuleRI::sternheimer_canonical_q_indices_one_based({q0, q1, q1_partner}),
              (std::vector<int>{1, 2}));

    q1.symmetry_spatial_isym = 2;
    EXPECT_THROW(ModuleRI::sternheimer_canonical_q_indices_one_based({q0, q1, q1_partner}),
                 std::invalid_argument);
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

TEST(SternheimerABACUSSTSmoke, BuildsDisjointFixedQKOrbitsFromLittleGroupPermutations)
{
    const std::vector<int> identity = {0, 1, 2, 3, 4, 5, 6, 7};
    const std::vector<int> inversion = {0, 7, 6, 5, 4, 3, 2, 1};

    const auto orbits
        = ModuleRI::build_sternheimer_fixed_q_k_orbits_from_permutations(8, {identity, inversion});

    ASSERT_EQ(orbits.size(), 5U);
    EXPECT_EQ(orbits[0].representative_ik_full, 0);
    EXPECT_EQ(orbits[0].members, (std::vector<int>{0}));
    EXPECT_EQ(orbits[1].representative_ik_full, 1);
    EXPECT_EQ(orbits[1].members, (std::vector<int>{1, 7}));
    EXPECT_EQ(orbits[2].members, (std::vector<int>{2, 6}));
    EXPECT_EQ(orbits[3].members, (std::vector<int>{3, 5}));
    EXPECT_EQ(orbits[4].members, (std::vector<int>{4}));
}

TEST(SternheimerABACUSSTSmoke, RejectsNonBijectiveFixedQOperation)
{
    EXPECT_THROW(ModuleRI::build_sternheimer_fixed_q_k_orbits_from_permutations(
                     4, {{0, 1, 2, 3}, {0, 0, 2, 3}}),
                 std::invalid_argument);
}

TEST(SternheimerABACUSSTSmoke, FormatsExplicitInverseRoutesForLibRPA)
{
    const std::vector<ModuleRI::SternheimerFixedQKRoute> routes = {
        {1, 3, 6, 2, true, {0, -1, 0}},
        {1, 0, 0, 0, false, {0, 0, 0}},
    };

    EXPECT_EQ(ModuleRI::format_sternheimer_fixed_q_routes(routes),
              "version 1\n"
              "# iq representative_ik member_ik spatial_isym time_reversal fold_Gx fold_Gy fold_Gz\n"
              "1 0 0 0 0 0 0 0\n"
              "1 3 6 2 1 0 -1 0\n");
    EXPECT_THROW(ModuleRI::format_sternheimer_fixed_q_routes({routes[0], routes[0]}),
                 std::invalid_argument);
}

TEST(SternheimerABACUSSTSmoke, BuildsDiscreteQStarRoutesFromAllowedGridPermutations)
{
    ModuleRI::SternheimerQStarPermutation identity;
    identity.spatial_isym = 0;
    identity.time_reversal = false;
    identity.mapped_index_by_full_q = {0, 1, 2, 3};
    identity.fold_G_by_full_q = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

    ModuleRI::SternheimerQStarPermutation swap;
    swap.spatial_isym = 4;
    swap.time_reversal = false;
    swap.mapped_index_by_full_q = {0, 2, 1, 3};
    swap.fold_G_by_full_q = {{0, 0, 0}, {-1, 0, 0}, {1, 0, 0}, {0, 0, 0}};

    const auto routes = ModuleRI::build_sternheimer_qstar_routes_from_permutations(
        4, {identity, swap});

    ASSERT_EQ(routes.size(), 4U);
    EXPECT_EQ(routes[0].representative_iq, 1);
    EXPECT_EQ(routes[0].member_iq, 1);
    EXPECT_EQ(routes[1].representative_iq, 2);
    EXPECT_EQ(routes[1].member_iq, 2);
    EXPECT_EQ(routes[2].representative_iq, 2);
    EXPECT_EQ(routes[2].member_iq, 3);
    EXPECT_EQ(routes[2].spatial_isym, 4);
    EXPECT_EQ(routes[2].fold_G, (std::array<int, 3>{1, 0, 0}));
    EXPECT_EQ(routes[3].representative_iq, 4);
    EXPECT_EQ(routes[3].member_iq, 4);

    EXPECT_EQ(ModuleRI::format_sternheimer_qstar_routes(routes),
              "version 1\n"
              "# representative_iq member_iq spatial_isym time_reversal fold_Gx fold_Gy fold_Gz\n"
              "1 1 0 0 0 0 0\n"
              "2 2 0 0 0 0 0\n"
              "2 3 4 0 1 0 0\n"
              "4 4 0 0 0 0 0\n");
}

TEST(SternheimerABACUSSTSmoke, HermitianizesEachPartialResponseBeforeOutput)
{
    using Complex = std::complex<double>;
    const std::vector<Complex> branch{{-1.0, 0.2}, {2.0, 3.0}, {4.0, -5.0}, {-6.0, 0.5}};

    const auto record = ModuleRI::make_sternheimer_partial_response_record(2, 3, 4, branch, 2);

    EXPECT_EQ(record.iq, 2);
    EXPECT_EQ(record.ik_full, 3);
    EXPECT_EQ(record.ifrequency, 4);
    EXPECT_EQ(record.filename, "v1_sternheimer_chi0_iq_2_ik_3_ifreq_4.dat");
    EXPECT_EQ(record.matrix,
              (std::vector<Complex>{{-2.0, 0.0}, {6.0, 8.0}, {6.0, -8.0}, {-12.0, 0.0}}));
}

TEST(SternheimerABACUSSTSmoke, SumOfPartialResponsesEqualsLegacyAggregateResponse)
{
    using Complex = std::complex<double>;
    const std::vector<Complex> branch_a{{-1.0, 0.2}, {2.0, 3.0}, {4.0, -5.0}, {-6.0, 0.5}};
    const std::vector<Complex> branch_b{{0.5, -0.2}, {-1.0, 1.0}, {2.0, 0.5}, {3.0, -0.5}};
    const auto partial_a
        = ModuleRI::make_sternheimer_partial_response_record(1, 0, 1, branch_a, 2);
    const auto partial_b
        = ModuleRI::make_sternheimer_partial_response_record(1, 1, 1, branch_b, 2);

    std::vector<Complex> partial_sum(partial_a.matrix.size(), Complex(0.0, 0.0));
    std::vector<Complex> aggregate_branch(branch_a.size(), Complex(0.0, 0.0));
    for (std::size_t index = 0; index != branch_a.size(); ++index)
    {
        partial_sum[index] = partial_a.matrix[index] + partial_b.matrix[index];
        aggregate_branch[index] = branch_a[index] + branch_b[index];
    }
    std::vector<Complex> expected(aggregate_branch.size(), Complex(0.0, 0.0));
    for (int row = 0; row != 2; ++row)
    {
        for (int column = 0; column != 2; ++column)
        {
            const std::size_t index = static_cast<std::size_t>(row * 2 + column);
            const std::size_t transpose = static_cast<std::size_t>(column * 2 + row);
            expected[index] = aggregate_branch[index] + std::conj(aggregate_branch[transpose]);
        }
    }
    EXPECT_EQ(partial_sum, expected);
}

TEST(SternheimerABACUSSTSmoke, FormatsPartialManifestInDeterministicKeyOrder)
{
    ModuleRI::SternheimerPartialResponseRecord later;
    later.iq = 2;
    later.ik_full = 7;
    later.ifrequency = 3;
    later.filename = "later.dat";
    ModuleRI::SternheimerPartialResponseRecord earlier;
    earlier.iq = 1;
    earlier.ik_full = 0;
    earlier.ifrequency = 1;
    earlier.filename = "earlier.dat";

    EXPECT_EQ(ModuleRI::format_sternheimer_partial_manifest({later, earlier}),
              "# iq ik_full ifreq response_file\n"
              "1 0 1 earlier.dat\n"
              "2 7 3 later.dat\n");
    EXPECT_THROW(ModuleRI::format_sternheimer_partial_manifest({earlier, earlier}),
                 std::invalid_argument);
}

TEST(SternheimerABACUSSTSmoke, FormatsFullKPointManifestInGlobalIndexOrder)
{
    const auto k2 = make_occupied_kpoint(2, 2, 0, {0.0, 0.5, 0.0}, 0.25);
    const auto k0 = make_occupied_kpoint(0, 0, 0, {0.0, 0.0, 0.0}, 0.25);
    const auto k1 = make_occupied_kpoint(1, 1, 0, {0.5, 0.0, 0.0}, 0.25);

    EXPECT_EQ(ModuleRI::format_sternheimer_full_kpoint_manifest({k2, k0, k1}),
              "# ik_full kx ky kz\n"
              "0 0 0 0\n"
              "1 0.5 0 0\n"
              "2 0 0.5 0\n");
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

TEST(SternheimerABACUSSTSmoke, BuildsSingleGammaSupercellTranslationPlanWithPositiveOutputIndex)
{
    const auto gamma = make_occupied_kpoint(0, 0, 0, {0.0, 0.0, 0.0}, 2.0);
    const auto plan = ModuleRI::build_sternheimer_periodic_response_plan({gamma}, 1, true);

    EXPECT_EQ(plan.iq, 1);
    EXPECT_EQ(plan.qpoint, (ModuleRI::SternheimerReducedKPoint{0.0, 0.0, 0.0}));
    ASSERT_EQ(plan.kq_pairs.size(), 1);
    EXPECT_EQ(plan.kq_pairs[0].source_index, 0);
    EXPECT_EQ(plan.kq_pairs[0].target_index, 0);
    EXPECT_DOUBLE_EQ(plan.kweight_sum, 2.0);
}

TEST(SternheimerABACUSSTSmoke, LimitsPeriodicOccupiedBandsOnlyWhenRequested)
{
    EXPECT_EQ(ModuleRI::sternheimer_periodic_band_count(8, -1), 8);
    EXPECT_EQ(ModuleRI::sternheimer_periodic_band_count(8, 1), 1);
    EXPECT_EQ(ModuleRI::sternheimer_periodic_band_count(8, 20), 8);
    EXPECT_THROW(ModuleRI::sternheimer_periodic_band_count(0, 1), std::invalid_argument);
}

TEST(SternheimerABACUSSTSmoke, TreatsSupercellTranslationResponseAsDiagnosticOnly)
{
    EXPECT_TRUE(ModuleRI::sternheimer_write_periodic_v1(false, false));
    EXPECT_FALSE(ModuleRI::sternheimer_write_periodic_v1(true, false));
    EXPECT_FALSE(ModuleRI::sternheimer_write_periodic_v1(false, true));

    EXPECT_NO_THROW(ModuleRI::validate_sternheimer_periodic_output_mode(true, true));
    EXPECT_NO_THROW(ModuleRI::validate_sternheimer_periodic_output_mode(false, false));
    EXPECT_THROW(ModuleRI::validate_sternheimer_periodic_output_mode(false, true),
                 std::invalid_argument);
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

TEST(SternheimerABACUSSTSmoke, AssignsEachNestedKFrequencyTaskExactlyOnce)
{
    constexpr int global_kpoint_count = 5;
    constexpr int kpoint_groups = 2;
    constexpr int frequency_count = 2;
    constexpr int mpi_ranks = kpoint_groups * frequency_count;
    std::vector<int> tasks_by_rank(mpi_ranks, 0);

    for (int ik = 0; ik != global_kpoint_count; ++ik)
    {
        for (int ifrequency = 0; ifrequency != frequency_count; ++ifrequency)
        {
            const auto assignment = ModuleRI::sternheimer_nested_mpi_assignment(ik,
                                                                                global_kpoint_count,
                                                                                ifrequency,
                                                                                frequency_count,
                                                                                kpoint_groups,
                                                                                mpi_ranks,
                                                                                0);
            EXPECT_EQ(assignment.kpoint_group,
                      ModuleRI::sternheimer_kpoint_owner_group(ik,
                                                               global_kpoint_count,
                                                               kpoint_groups));
            EXPECT_EQ(assignment.frequency_slot, ifrequency);
            ASSERT_GE(assignment.owner_rank, 0);
            ASSERT_LT(assignment.owner_rank, mpi_ranks);
            ++tasks_by_rank[static_cast<std::size_t>(assignment.owner_rank)];
        }
    }

    EXPECT_EQ(tasks_by_rank, (std::vector<int>{3, 3, 2, 2}));
}

TEST(SternheimerABACUSSTSmoke, WrapsNestedFrequencyRankShiftWithinEachKGroup)
{
    const auto positive = ModuleRI::sternheimer_nested_mpi_assignment(3, 5, 1, 3, 2, 6, 1);
    EXPECT_EQ(positive.kpoint_group, 1);
    EXPECT_EQ(positive.frequency_slot, 2);
    EXPECT_EQ(positive.owner_rank, 5);

    const auto negative = ModuleRI::sternheimer_nested_mpi_assignment(0, 5, 0, 3, 2, 6, -1);
    EXPECT_EQ(negative.kpoint_group, 0);
    EXPECT_EQ(negative.frequency_slot, 2);
    EXPECT_EQ(negative.owner_rank, 2);
}

TEST(SternheimerABACUSSTSmoke, RejectsInvalidNestedMPIContracts)
{
    EXPECT_THROW(ModuleRI::sternheimer_nested_mpi_assignment(0, 4, 0, 2, 2, 3),
                 std::invalid_argument);
    EXPECT_THROW(ModuleRI::sternheimer_nested_mpi_assignment(-1, 4, 0, 2, 2, 4),
                 std::invalid_argument);
    EXPECT_THROW(ModuleRI::sternheimer_nested_mpi_assignment(0, 4, -1, 2, 2, 4),
                 std::invalid_argument);
    EXPECT_THROW(ModuleRI::sternheimer_nested_mpi_assignment(0, 4, 2, 2, 2, 4),
                 std::invalid_argument);
    EXPECT_THROW(ModuleRI::sternheimer_nested_mpi_assignment(0, 2, 0, 2, 3, 6),
                 std::invalid_argument);
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
