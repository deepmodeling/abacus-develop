#include "source_lcao/module_ri/sternheimer_abacus_st_smoke.h"

#include <complex>
#include <gtest/gtest.h>

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

TEST(SternheimerABACUSSTSmoke, ValidatesSpinResolvedLCAOOccupiedChannels)
{
    ModuleRI::SternheimerLCAOOccupiedChannel spin_up;
    spin_up.spin_index = 0;
    spin_up.coefficients = {{std::complex<double>(1.0, 0.0),
                             std::complex<double>(0.0, 0.0),
                             std::complex<double>(0.0, 0.0)}};

    const std::vector<ModuleRI::SternheimerLCAOOccupiedChannel> channels = {spin_up};
    EXPECT_NO_THROW(ModuleRI::validate_sternheimer_lcao_occupied_channels(channels, 2, 3));
    EXPECT_EQ(ModuleRI::sternheimer_lcao_total_occupied_bands(channels), 1);

    auto duplicate_channels = channels;
    duplicate_channels.push_back(spin_up);
    EXPECT_THROW(ModuleRI::validate_sternheimer_lcao_occupied_channels(duplicate_channels, 2, 3),
                 std::invalid_argument);

    auto out_of_range_channels = channels;
    out_of_range_channels.front().spin_index = 2;
    EXPECT_THROW(ModuleRI::validate_sternheimer_lcao_occupied_channels(out_of_range_channels, 2, 3),
                 std::invalid_argument);
}
