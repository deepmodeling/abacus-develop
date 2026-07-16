#include "source_lcao/module_ri/sternheimer_abacus_st_smoke.h"

#include <complex>
#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>

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

TEST(SternheimerABACUSSTSmoke, RejectsGammaOrFractionalOccupationForSolidQ)
{
    auto k0 = make_occupied_kpoint(0, 0, 0, {0.0, 0.0, 0.0}, 1.0);
    const auto k1 = make_occupied_kpoint(1, 1, 0, {0.5, 0.0, 0.0}, 1.0);

    EXPECT_THROW(ModuleRI::build_sternheimer_periodic_response_plan({k0, k1}, 1), std::invalid_argument);
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
