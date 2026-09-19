#include "source_lcao/module_ri/sternheimer_abacus_st_smoke.h"

#include "source_lcao/module_ri/sternheimer_channel_parallel.h"

#include <array>
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

TEST(SternheimerABACUSSTSmoke, UsesProductionDefaults)
{
    EXPECT_DOUBLE_EQ(ModuleRI::default_sternheimer_solver_tolerance(), 1.0e-6);
    EXPECT_EQ(ModuleRI::parse_sternheimer_lcao_virtual_source(""),
              ModuleRI::SternheimerLCAOVirtualSource::KSBands);
}

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

#ifdef _OPENMP
    omp_set_num_threads(previous_threads);
    omp_set_dynamic(previous_dynamic);
#endif
    ASSERT_EQ(results.size(), 32U);
    for (int channel_index = 0; channel_index != 32; ++channel_index)
    {
        EXPECT_EQ(results[static_cast<std::size_t>(channel_index)], channel_index * channel_index);
    }
#ifdef _OPENMP
    EXPECT_GT(peak_workers.load(), 1);
#else
    EXPECT_EQ(peak_workers.load(), 1);
#endif
}

TEST(SternheimerChannelParallel, HonorsExplicitMaximumWorkerCount)
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
        = ModuleRI::run_sternheimer_channel_tasks<int>(
            32,
            [&active_workers, &peak_workers](const int channel_index) {
                const int active = active_workers.fetch_add(1) + 1;
                int observed_peak = peak_workers.load();
                while (active > observed_peak && !peak_workers.compare_exchange_weak(observed_peak, active))
                {
                }
                std::this_thread::sleep_for(std::chrono::milliseconds(2));
                active_workers.fetch_sub(1);
                return channel_index;
            },
            2);

#ifdef _OPENMP
    omp_set_num_threads(previous_threads);
    omp_set_dynamic(previous_dynamic);
#endif
    ASSERT_EQ(results.size(), 32U);
    for (int channel_index = 0; channel_index != 32; ++channel_index)
    {
        EXPECT_EQ(results[static_cast<std::size_t>(channel_index)], channel_index);
    }
#ifdef _OPENMP
    EXPECT_EQ(peak_workers.load(), 2);
#else
    EXPECT_EQ(peak_workers.load(), 1);
#endif
}

TEST(SternheimerChannelParallel, SingleWorkerPreservesNestedGridTeam)
{
#ifdef _OPENMP
    const int previous_threads = omp_get_max_threads();
    const int previous_dynamic = omp_get_dynamic();
    const int previous_active_levels = omp_get_max_active_levels();
    omp_set_dynamic(0);
    omp_set_num_threads(4);
    omp_set_max_active_levels(1);
#endif

    const std::vector<int> nested_team_sizes = ModuleRI::run_sternheimer_channel_tasks<int>(
        4,
        [](const int) {
#ifdef _OPENMP
            int nested_team_size = 0;
#pragma omp parallel
            {
#pragma omp single
                nested_team_size = omp_get_num_threads();
            }
            return nested_team_size;
#else
            return 1;
#endif
        },
        1);

#ifdef _OPENMP
    omp_set_num_threads(previous_threads);
    omp_set_dynamic(previous_dynamic);
    omp_set_max_active_levels(previous_active_levels);
#endif
    ASSERT_EQ(nested_team_sizes.size(), 4U);
    for (const int nested_team_size: nested_team_sizes)
    {
#ifdef _OPENMP
        EXPECT_EQ(nested_team_size, 4);
#else
        EXPECT_EQ(nested_team_size, 1);
#endif
    }
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

    spin_up.coefficients.resize(4, spin_up.coefficients.front());
    ModuleRI::SternheimerLCAOOccupiedChannel spin_down;
    spin_down.spin_index = 1;
    spin_down.coefficients = {{std::complex<double>(0.0, 0.0),
                               std::complex<double>(1.0, 0.0),
                               std::complex<double>(0.0, 0.0)}};
    const std::vector<ModuleRI::SternheimerLCAOOccupiedChannel> quartet = {spin_up, spin_down};

    EXPECT_NO_THROW(ModuleRI::validate_sternheimer_lcao_occupied_channels(quartet, 2, 3));
    EXPECT_EQ(ModuleRI::sternheimer_lcao_spin_indices(quartet), (std::vector<int>{0, 1}));
    EXPECT_EQ(ModuleRI::sternheimer_lcao_occupied_bands_per_spin(quartet), (std::vector<int>{4, 1}));
    EXPECT_EQ(ModuleRI::sternheimer_lcao_total_occupied_bands(quartet), 5);

    auto duplicate_channels = channels;
    duplicate_channels.push_back(spin_up);
    EXPECT_THROW(ModuleRI::validate_sternheimer_lcao_occupied_channels(duplicate_channels, 2, 3),
                 std::invalid_argument);

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

TEST(SternheimerABACUSSTSmoke, UnlimitedDiagnosticChannelsReuseRvalueStorage)
{
    std::vector<ModuleRI::SternheimerABFBlochGridChannel> channels(1);
    channels[0].potential_r.resize(8, std::complex<double>(1.0, -0.5));
    const auto* original_storage = channels[0].potential_r.data();

    const auto unlimited = ModuleRI::limit_sternheimer_abf_channels_per_atom(std::move(channels), -1);

    ASSERT_EQ(unlimited.size(), 1);
    EXPECT_EQ(unlimited[0].potential_r.data(), original_storage);
}

TEST(SternheimerABACUSSTSmoke, SamplesOnlyRequestedPeriodicKSUnoccupiedStates)
{
    EXPECT_EQ(ModuleRI::sternheimer_sampled_unoccupied_count(false, 676, 32), 0);
    EXPECT_EQ(ModuleRI::sternheimer_sampled_unoccupied_count(true, 676, 0), 676);
    EXPECT_EQ(ModuleRI::sternheimer_sampled_unoccupied_count(true, 676, 32), 32);
    EXPECT_EQ(ModuleRI::sternheimer_sampled_unoccupied_count(true, 12, 32), 12);
    EXPECT_THROW(ModuleRI::sternheimer_sampled_unoccupied_count(true, 676, -1), std::invalid_argument);
}

TEST(SternheimerABACUSSTSmoke, SamplesOnlyRequestedPeriodicOccupiedSourceStates)
{
    EXPECT_EQ(ModuleRI::sternheimer_sampled_occupied_count(68, 1), 1);
    EXPECT_EQ(ModuleRI::sternheimer_sampled_occupied_count(68, 0), 68);
    EXPECT_EQ(ModuleRI::sternheimer_sampled_occupied_count(68, -1), 68);
    EXPECT_EQ(ModuleRI::sternheimer_sampled_occupied_count(12, 32), 12);
    EXPECT_THROW(ModuleRI::sternheimer_sampled_occupied_count(0, 1), std::invalid_argument);
}

TEST(SternheimerABACUSSTSmoke, GathersRequestedPeriodicKSUnoccupiedStatesForDelta)
{
    EXPECT_EQ(ModuleRI::sternheimer_lcao_virtual_state_gather_count(676, false, "ks_bands", 32), 0);
    EXPECT_EQ(ModuleRI::sternheimer_lcao_virtual_state_gather_count(676, true, "projected_ao", 32), 0);
    EXPECT_EQ(ModuleRI::sternheimer_lcao_virtual_state_gather_count(676, true, "ks_bands", 0), 676);
    EXPECT_EQ(ModuleRI::sternheimer_lcao_virtual_state_gather_count(676, true, "ks_bands", 32), 32);
    EXPECT_EQ(ModuleRI::sternheimer_lcao_virtual_state_gather_count(12, true, "ks_bands", 32), 12);
    EXPECT_THROW(ModuleRI::sternheimer_lcao_virtual_state_gather_count(-1, true, "ks_bands", 32),
                 std::invalid_argument);
    EXPECT_THROW(ModuleRI::sternheimer_lcao_virtual_state_gather_count(676, true, "ks_bands", -1),
                 std::invalid_argument);
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

TEST(SternheimerABACUSSTSmoke, DoesNotReconstructFullKGridWhenSymmetryIsDisabled)
{
    EXPECT_FALSE(ModuleRI::sternheimer_full_k_reconstruction_required(2, 1, 2, 0));
}

TEST(SternheimerABACUSSTSmoke, DoesNotReconstructAnAlreadyFullKGrid)
{
    EXPECT_FALSE(ModuleRI::sternheimer_full_k_reconstruction_required(4, 4, 1, 1));
}

TEST(SternheimerABACUSSTSmoke, DoesNotReconstructFullGammaForTwoSpinChannels)
{
    EXPECT_FALSE(ModuleRI::sternheimer_full_k_reconstruction_required(2, 1, 2, 1));
}

TEST(SternheimerABACUSSTSmoke, ReconstructsASymmetryReducedSingleSpinKGrid)
{
    EXPECT_TRUE(ModuleRI::sternheimer_full_k_reconstruction_required(2, 8, 1, 1));
}

TEST(SternheimerABACUSSTSmoke, RejectsSymmetryReducedSpinPolarizedKGrid)
{
    EXPECT_THROW(ModuleRI::sternheimer_full_k_reconstruction_required(4, 8, 2, 1),
                 std::invalid_argument);
}

TEST(SternheimerABACUSSTSmoke, SelectsExplicitKSVirtualBandsFromOccupations)
{
    ModuleRI::SternheimerLCAOOccupiedChannel spin_up;
    spin_up.spin_index = 0;
    spin_up.coefficients = {{std::complex<double>(1.0, 0.0),
                             std::complex<double>(0.0, 0.0),
                             std::complex<double>(0.0, 0.0)}};
    spin_up.unoccupied_coefficients = {
        {std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0)},
        {std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0)}};

    const std::vector<ModuleRI::SternheimerLCAOOccupiedChannel> channels = {spin_up};
    EXPECT_NO_THROW(ModuleRI::validate_sternheimer_lcao_occupied_channels(channels, 1, 3));
    EXPECT_EQ(ModuleRI::sternheimer_lcao_unoccupied_bands_per_spin(channels), (std::vector<int>{2}));
    EXPECT_EQ(ModuleRI::sternheimer_lcao_total_unoccupied_bands(channels), 2);

    EXPECT_EQ(ModuleRI::parse_sternheimer_lcao_virtual_source("projected_ao"),
              ModuleRI::SternheimerLCAOVirtualSource::ProjectedAO);
    EXPECT_EQ(ModuleRI::parse_sternheimer_lcao_virtual_source("ks_bands"),
              ModuleRI::SternheimerLCAOVirtualSource::KSBands);
    EXPECT_EQ(ModuleRI::sternheimer_lcao_virtual_source_name(ModuleRI::SternheimerLCAOVirtualSource::KSBands),
              "ks_bands");
    EXPECT_THROW(ModuleRI::parse_sternheimer_lcao_virtual_source("svd_guess"), std::invalid_argument);

    auto invalid = channels;
    invalid.front().unoccupied_coefficients.front().pop_back();
    EXPECT_THROW(ModuleRI::validate_sternheimer_lcao_occupied_channels(invalid, 1, 3), std::invalid_argument);
}

TEST(SternheimerABACUSSTSmoke, RequiresCompleteKSVirtualSubspace)
{
    EXPECT_EQ(ModuleRI::expected_sternheimer_ks_virtual_states(21, 0), 21);
    EXPECT_EQ(ModuleRI::expected_sternheimer_ks_virtual_states(21, 16), 16);
    EXPECT_NO_THROW(ModuleRI::validate_sternheimer_ks_virtual_subspace(1, 21, 21, 21));
    EXPECT_THROW(ModuleRI::validate_sternheimer_ks_virtual_subspace(2, 21, 20, 20),
                 std::runtime_error);
    EXPECT_THROW(ModuleRI::validate_sternheimer_ks_virtual_subspace(2, 21, 21, 20),
                 std::runtime_error);
}

TEST(SternheimerABACUSSTSmoke, AcceptsOnlyPhysicalGammaSpinRows)
{
    const std::vector<std::array<double, 3>> one_gamma = {{{0.0, 0.0, 0.0}}};
    const std::vector<std::array<double, 3>> two_gamma = {{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}};

    EXPECT_NO_THROW(ModuleRI::validate_sternheimer_lcao_gamma_layout(1, 1, 1, one_gamma));
    EXPECT_NO_THROW(ModuleRI::validate_sternheimer_lcao_gamma_layout(2, 2, 2, two_gamma));

    EXPECT_THROW(ModuleRI::validate_sternheimer_lcao_gamma_layout(1, 2, 2, two_gamma), std::invalid_argument);

    const std::vector<std::array<double, 3>> non_gamma = {{{0.25, 0.0, 0.0}}};
    EXPECT_THROW(ModuleRI::validate_sternheimer_lcao_gamma_layout(1, 1, 1, non_gamma), std::invalid_argument);
}

TEST(SternheimerABACUSSTSmoke, SelectsZeroOrderSourceByResponseMode)
{
    EXPECT_FALSE(ModuleRI::sternheimer_uses_lcao_zero_order(false));
    EXPECT_TRUE(ModuleRI::sternheimer_uses_lcao_zero_order(true));
}

TEST(SternheimerABACUSSTSmoke, MapsCanonicalFullQIndexToReaderV1StarIndex)
{
    auto q0 = make_occupied_kpoint(0, 0, 0, {0.0, 0.0, 0.0}, 0.25);
    auto qx = make_occupied_kpoint(1, 1, 0, {0.25, 0.0, 0.0}, 0.25);
    auto qx_partner = make_occupied_kpoint(2, 2, 0, {-0.25, 0.0, 0.0}, 0.25);
    auto qy = make_occupied_kpoint(3, 3, 0, {0.0, 0.25, 0.0}, 0.25);
    qx_partner.zero_order_k_index = 1;
    qx_partner.symmetry_spatial_isym = 2;
    qy.zero_order_k_index = 2;

    const std::vector<ModuleRI::SternheimerLCAOOccupiedKPoint> records{
        q0, qx, qx_partner, qy};
    EXPECT_EQ(ModuleRI::sternheimer_coulomb_reader_q_index_one_based(records, 1), 1);
    EXPECT_EQ(ModuleRI::sternheimer_coulomb_reader_q_index_one_based(records, 2), 2);
    EXPECT_EQ(ModuleRI::sternheimer_coulomb_reader_q_index_one_based(records, 4), 3);
    EXPECT_THROW(ModuleRI::sternheimer_coulomb_reader_q_index_one_based(records, 0),
                 std::invalid_argument);
}

TEST(SternheimerABACUSSTSmoke, SelectsEveryGammaSpinRecordForMolecularResponse)
{
    EXPECT_EQ(ModuleRI::sternheimer_abfs_perturbation_source({}), "product_pca");
    EXPECT_EQ(ModuleRI::sternheimer_abfs_perturbation_source({"H-fixed.abfs"}), "explicit_abfs");
    EXPECT_TRUE(ModuleRI::sternheimer_builds_product_pca_auxiliary_basis({}));
    EXPECT_FALSE(ModuleRI::sternheimer_builds_product_pca_auxiliary_basis({"H-fixed.abfs"}));
}

TEST(SternheimerABACUSSTSmoke, SelectsOnlyOccupiedGammaSpinRecordsForOpenShellResponse)
{
    auto occupied_spin = make_occupied_kpoint(1, 1, 1, {0.0, 0.0, 0.0}, 1.0);

    EXPECT_NO_THROW(ModuleRI::validate_sternheimer_lcao_occupied_kpoints(
        {occupied_spin}, 2, 2, 2, 3, -1, false));
    const auto selected = ModuleRI::select_sternheimer_gamma_spin_records({occupied_spin}, 2);

    ASSERT_EQ(selected.size(), 1U);
    EXPECT_EQ(selected.front()->spin_index, 1);
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

TEST(SternheimerABACUSSTSmoke, DividesSupercellKWeightAcrossPrimitiveTranslationSectors)
{
    EXPECT_DOUBLE_EQ(ModuleRI::sternheimer_supercell_sector_kweight(2.0, 64), 0.03125);
    EXPECT_DOUBLE_EQ(ModuleRI::sternheimer_supercell_sector_kweight(0.5, 8), 1.0 / 16.0);
    EXPECT_THROW(ModuleRI::sternheimer_supercell_sector_kweight(1.0, 0), std::invalid_argument);
}

TEST(SternheimerABACUSSTSmoke, RestoresPrimitiveNormalizationForFullSupercellResponseMatrix)
{
    EXPECT_DOUBLE_EQ(ModuleRI::sternheimer_supercell_response_matrix_scale(false, 64), 1.0);
    EXPECT_DOUBLE_EQ(ModuleRI::sternheimer_supercell_response_matrix_scale(true, 8), 8.0);
    EXPECT_DOUBLE_EQ(ModuleRI::sternheimer_supercell_response_matrix_scale(true, 64), 64.0);
    EXPECT_THROW(ModuleRI::sternheimer_supercell_response_matrix_scale(true, 0),
                 std::invalid_argument);
}

TEST(SternheimerABACUSSTSmoke, PlansVirtualSamplingWithoutMakingTheOptionalSOSDiagnosticMandatory)
{
    const auto standard = ModuleRI::sternheimer_lcao_sampling_plan(false, false, false);
    EXPECT_FALSE(standard.sample_source_unoccupied);
    EXPECT_FALSE(standard.sample_target_unoccupied);
    EXPECT_FALSE(standard.build_target_ao_candidates);

    const auto delta_fallback = ModuleRI::sternheimer_lcao_sampling_plan(true, false, false);
    EXPECT_FALSE(delta_fallback.sample_source_unoccupied);
    EXPECT_FALSE(delta_fallback.sample_target_unoccupied);
    EXPECT_TRUE(delta_fallback.build_target_ao_candidates);

    const auto delta_virtuals = ModuleRI::sternheimer_lcao_sampling_plan(true, false, true);
    EXPECT_FALSE(delta_virtuals.sample_source_unoccupied);
    EXPECT_TRUE(delta_virtuals.sample_target_unoccupied);
    EXPECT_FALSE(delta_virtuals.build_target_ao_candidates);

    EXPECT_THROW(ModuleRI::sternheimer_lcao_sampling_plan(false, true, false), std::invalid_argument);
}

TEST(SternheimerABACUSSTSmoke, ReusesTargetSamplingWhenItContainsTheSameSourceRecord)
{
    const auto occupied_only = ModuleRI::sternheimer_lcao_sampling_plan(false, false, false);
    const auto target_virtuals = ModuleRI::sternheimer_lcao_sampling_plan(true, false, true);
    ModuleRI::SternheimerLCAOSamplingPlan missing_source_virtuals;
    missing_source_virtuals.sample_source_unoccupied = true;

    EXPECT_TRUE(ModuleRI::can_reuse_sternheimer_target_lcao_sampling(true, occupied_only));
    EXPECT_TRUE(ModuleRI::can_reuse_sternheimer_target_lcao_sampling(true, target_virtuals));
    EXPECT_FALSE(ModuleRI::can_reuse_sternheimer_target_lcao_sampling(false, target_virtuals));
    EXPECT_FALSE(ModuleRI::can_reuse_sternheimer_target_lcao_sampling(true, missing_source_virtuals));
}

TEST(SternheimerABACUSSTSmoke, RejectsNonInsulatingSupercellSectorOccupations)
{
    EXPECT_NO_THROW(ModuleRI::validate_sternheimer_supercell_sector_occupations({1.0, 1.0}, 2));
    EXPECT_THROW(ModuleRI::validate_sternheimer_supercell_sector_occupations({1.0}, 2),
                 std::invalid_argument);
    EXPECT_THROW(ModuleRI::validate_sternheimer_supercell_sector_occupations({1.0, 0.5}, 2),
                 std::invalid_argument);
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

TEST(SternheimerABACUSSTSmoke, BuildsSingleGammaPeriodicResponsePlanWithPositiveIndex)
{
    const auto gamma = make_occupied_kpoint(0, 0, 0, {0.0, 0.0, 0.0}, 2.0);
    const auto plan = ModuleRI::build_sternheimer_periodic_response_plan({gamma}, 1);

    EXPECT_EQ(plan.iq, 1);
    EXPECT_EQ(plan.qpoint, (ModuleRI::SternheimerReducedKPoint{0.0, 0.0, 0.0}));
    EXPECT_EQ(plan.record_index_by_global_k, (std::vector<int>{0}));
    ASSERT_EQ(plan.kq_pairs.size(), 1U);
    EXPECT_EQ(plan.kq_pairs[0].source_index, 0);
    EXPECT_EQ(plan.kq_pairs[0].target_index, 0);
    EXPECT_DOUBLE_EQ(plan.kweight_sum, 2.0);
    EXPECT_NO_THROW(ModuleRI::validate_sternheimer_periodic_kmesh({1, 1, 1}, 1));
}

TEST(SternheimerABACUSSTSmoke, RejectsNonGammaSinglePointPeriodicResponse)
{
    for (int direction = 0; direction != 3; ++direction)
    {
        auto point = make_occupied_kpoint(0, 0, 0, {0.0, 0.0, 0.0}, 2.0);
        point.kpoint[direction] = 0.25;
        EXPECT_THROW(ModuleRI::build_sternheimer_periodic_response_plan({point}, 1),
                     std::invalid_argument);
    }
}

TEST(SternheimerABACUSSTSmoke, RejectsFractionalOccupationForSingleGammaPeriodicResponse)
{
    auto gamma = make_occupied_kpoint(0, 0, 0, {0.0, 0.0, 0.0}, 2.0);
    gamma.occupations[0] = 0.5;
    EXPECT_THROW(ModuleRI::build_sternheimer_periodic_response_plan({gamma}, 1),
                 std::invalid_argument);
}

TEST(SternheimerABACUSSTSmoke, RejectsSpinChannelForSingleGammaPeriodicResponse)
{
    const auto gamma = make_occupied_kpoint(0, 0, 1, {0.0, 0.0, 0.0}, 1.0);
    EXPECT_THROW(ModuleRI::build_sternheimer_periodic_response_plan({gamma}, 1),
                 std::invalid_argument);
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

TEST(SternheimerABACUSSTSmoke, SelectsWeakQRepresentativeAndCompleteInverseRoutes)
{
    const std::vector<ModuleRI::SternheimerFixedQKOrbit> orbits = {
        {0, {0}},
        {1, {1, 3}},
        {2, {2}},
    };
    const std::vector<ModuleRI::SternheimerFixedQKRoute> routes = {
        {2, 0, 0, 0, false, {0, 0, 0}},
        {2, 1, 1, 0, false, {0, 0, 0}},
        {2, 1, 3, 4, true, {1, 0, 0}},
        {2, 2, 2, 0, false, {0, 0, 0}},
    };

    const auto representative
        = ModuleRI::select_sternheimer_weak_q_symmetry_source(2, 2, orbits, routes);
    EXPECT_EQ(representative.iq, 2);
    EXPECT_EQ(representative.source_ik_full, 1);
    EXPECT_EQ(representative.representative_ik_full, 1);
    EXPECT_TRUE(representative.source_is_representative);
    EXPECT_EQ(representative.orbit_size, 2);
    ASSERT_EQ(representative.inverse_routes.size(), 2U);
    EXPECT_EQ(representative.inverse_routes[0].member_ik_full, 1);
    EXPECT_EQ(representative.inverse_routes[1].member_ik_full, 3);

    const auto member
        = ModuleRI::select_sternheimer_weak_q_symmetry_source(2, 4, orbits, routes);
    EXPECT_EQ(member.source_ik_full, 3);
    EXPECT_EQ(member.representative_ik_full, 1);
    EXPECT_FALSE(member.source_is_representative);
    EXPECT_EQ(member.orbit_size, 2);
}

TEST(SternheimerABACUSSTSmoke, RejectsIncompleteWeakQSymmetryRoutes)
{
    const std::vector<ModuleRI::SternheimerFixedQKOrbit> orbits = {
        {0, {0, 2}},
        {1, {1}},
    };
    const std::vector<ModuleRI::SternheimerFixedQKRoute> missing_member = {
        {1, 0, 0, 0, false, {0, 0, 0}},
        {1, 1, 1, 0, false, {0, 0, 0}},
    };
    EXPECT_THROW(ModuleRI::select_sternheimer_weak_q_symmetry_source(
                     1, 1, orbits, missing_member),
                 std::invalid_argument);

    auto duplicate_member = missing_member;
    duplicate_member.push_back({1, 0, 0, 0, false, {0, 0, 0}});
    duplicate_member.push_back({1, 0, 2, 3, false, {0, 0, 0}});
    EXPECT_THROW(ModuleRI::select_sternheimer_weak_q_symmetry_source(
                     1, 1, orbits, duplicate_member),
                 std::invalid_argument);
}

TEST(SternheimerABACUSSTSmoke, WeakQSymmetryRestorationRequiresAllOccupiedBands)
{
    EXPECT_NO_THROW(ModuleRI::validate_sternheimer_weak_q_symmetry_band_coverage(
        1, 5, 5));
    EXPECT_THROW(ModuleRI::validate_sternheimer_weak_q_symmetry_band_coverage(
                     1, 3, 5),
                 std::invalid_argument);
    EXPECT_THROW(ModuleRI::validate_sternheimer_weak_q_symmetry_band_coverage(
                     2, 5, 5),
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

TEST(SternheimerABACUSSTSmoke, ExpandsGammaSupercellIntoPrimitiveKRecordsWithGammaGridBoundary)
{
    ModuleRI::SternheimerLCAOOccupiedKPoint gamma;
    gamma.local_k_index = 0;
    gamma.global_k_index = 0;
    gamma.zero_order_k_index = 0;
    gamma.spin_index = 0;
    gamma.kpoint = {0.0, 0.0, 0.0};
    gamma.kweight = 2.0;
    gamma.eigenvalues = {-2.0, -1.0};
    gamma.occupations = {1.0, 1.0};
    gamma.coefficients = {
        {1.0, 0.0, 0.0, 0.0},
        {0.0, 1.0, 0.0, 0.0},
    };
    gamma.unoccupied_eigenvalues = {1.0, 2.0};
    gamma.unoccupied_coefficients = {
        {0.0, 0.0, 1.0, 0.0},
        {0.0, 0.0, 0.0, 1.0},
    };

    std::vector<ModuleRI::SternheimerSupercellKPointSector> sectors(2);
    sectors[0].kpoint = {0.0, 0.0, 0.0};
    sectors[0].sector.eigenvalues = {-2.0, 1.0};
    sectors[0].sector.coefficients = {gamma.coefficients[0], gamma.unoccupied_coefficients[0]};
    sectors[1].kpoint = {0.5, 0.0, 0.0};
    sectors[1].sector.eigenvalues = {-1.0, 2.0};
    sectors[1].sector.coefficients = {gamma.coefficients[1], gamma.unoccupied_coefficients[1]};

    const auto records
        = ModuleRI::build_sternheimer_supercell_full_kpoint_records(gamma, sectors);

    ASSERT_EQ(records.size(), 2U);
    for (std::size_t ik = 0; ik != records.size(); ++ik)
    {
        EXPECT_EQ(records[ik].local_k_index, static_cast<int>(ik));
        EXPECT_EQ(records[ik].global_k_index, static_cast<int>(ik));
        EXPECT_EQ(records[ik].zero_order_k_index, static_cast<int>(ik));
        EXPECT_DOUBLE_EQ(records[ik].kweight, 1.0);
        EXPECT_EQ(records[ik].kpoint, sectors[ik].kpoint);
        EXPECT_EQ(ModuleRI::sternheimer_lcao_grid_kpoint(records[ik]),
                  (ModuleRI::SternheimerReducedKPoint{0.0, 0.0, 0.0}));
        ASSERT_EQ(records[ik].eigenvalues.size(), 1U);
        ASSERT_EQ(records[ik].unoccupied_eigenvalues.size(), 1U);
        EXPECT_DOUBLE_EQ(records[ik].occupations[0], 1.0);
    }
    EXPECT_DOUBLE_EQ(records[0].eigenvalues[0], -2.0);
    EXPECT_DOUBLE_EQ(records[1].eigenvalues[0], -1.0);
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
    EXPECT_TRUE(ModuleRI::sternheimer_write_periodic_v1(false, true, false, true));
    EXPECT_TRUE(ModuleRI::sternheimer_write_periodic_v1(true, false, true));
    EXPECT_FALSE(ModuleRI::sternheimer_write_periodic_v1(true, true, true));
    EXPECT_TRUE(ModuleRI::sternheimer_write_periodic_v1(true, true, true, true));

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

TEST(SternheimerABACUSSTSmoke, UsesTheTwoDimensionalMassiddaIntegralForTwoDimensionalGamma)
{
    EXPECT_TRUE(ModuleRI::sternheimer_periodic_gamma_uses_2d_massidda({0.0, 0.0, 0.0}, 2));
    EXPECT_FALSE(ModuleRI::sternheimer_periodic_gamma_uses_2d_massidda({0.25, 0.0, 0.0}, 2));
    EXPECT_FALSE(ModuleRI::sternheimer_periodic_gamma_uses_2d_massidda({0.0, 0.0, 0.0}, 3));
    EXPECT_THROW(ModuleRI::sternheimer_periodic_gamma_uses_2d_massidda({0.0, 0.0, 0.0}, 1),
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

TEST(SternheimerABACUSSTSmoke, SelectsIndependentSupercellKPointGroups)
{
    EXPECT_EQ(ModuleRI::sternheimer_response_kpoint_group_count(false, 1, 4, 64), 4);
    EXPECT_EQ(ModuleRI::sternheimer_response_kpoint_group_count(true, 8, 1, 64), 8);
    EXPECT_EQ(ModuleRI::sternheimer_response_kpoint_group_count(true, 64, 1, 64), 64);
    EXPECT_THROW(ModuleRI::sternheimer_response_kpoint_group_count(true, 0, 1, 64),
                 std::invalid_argument);
    EXPECT_THROW(ModuleRI::sternheimer_response_kpoint_group_count(true, 65, 1, 64),
                 std::invalid_argument);
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

TEST(SternheimerABACUSSTSmoke, BuildsNestedChannelMPIReplicaLayout)
{
    const auto base = ModuleRI::sternheimer_nested_mpi_replica_layout(2, 16, 32, 31, false);
    EXPECT_EQ(base.response_slots, 32);
    EXPECT_EQ(base.replicas_per_slot, 1);
    EXPECT_EQ(base.local_response_slot, 31);
    EXPECT_EQ(base.local_replica, 0);

    const auto replicated = ModuleRI::sternheimer_nested_mpi_replica_layout(2, 16, 128, 127, true);
    EXPECT_EQ(replicated.response_slots, 32);
    EXPECT_EQ(replicated.replicas_per_slot, 4);
    EXPECT_EQ(replicated.local_response_slot, 31);
    EXPECT_EQ(replicated.local_replica, 3);

    EXPECT_THROW(ModuleRI::sternheimer_nested_mpi_replica_layout(2, 16, 64, 0, false),
                 std::invalid_argument);
    EXPECT_THROW(ModuleRI::sternheimer_nested_mpi_replica_layout(2, 16, 48, 0, true),
                 std::invalid_argument);
}

TEST(SternheimerABACUSSTSmoke, AssignsChannelBatchesOnceAndBalancesAcrossReplicas)
{
    constexpr int occupied_states = 5;
    constexpr int batch_count = 126;
    constexpr int replica_count = 8;
    std::vector<int> tasks_by_replica(replica_count, 0);
    for (int occupied = 0; occupied != occupied_states; ++occupied)
    {
        for (int batch = 0; batch != batch_count; ++batch)
        {
            const int owner = ModuleRI::sternheimer_channel_batch_replica_owner(
                occupied, batch, batch_count, replica_count);
            ASSERT_GE(owner, 0);
            ASSERT_LT(owner, replica_count);
            ++tasks_by_replica[static_cast<std::size_t>(owner)];
        }
    }
    const auto range = std::minmax_element(tasks_by_replica.begin(), tasks_by_replica.end());
    EXPECT_LE(*range.second - *range.first, 1);
    EXPECT_THROW(ModuleRI::sternheimer_channel_batch_replica_owner(0, 126, 126, 8),
                 std::invalid_argument);
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
