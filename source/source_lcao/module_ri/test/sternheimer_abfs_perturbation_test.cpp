#include "source_lcao/module_ri/sternheimer_abfs_perturbation.h"
#include "source_lcao/module_ri/sternheimer_abacus_st_smoke.h"

#include <gtest/gtest.h>

#include <cstdlib>

TEST(SternheimerABFSPerturbation, SamplesSChannelOnFDGrid)
{
    ModuleRI::SternheimerRadialPerturbation radial;
    radial.angular_momentum = 0;
    radial.radial_index = 0;
    radial.label = "toy_s";
    radial.radial_grid = {0.0, 1.0};
    radial.radial_values = {2.0, 4.0};

    ModuleRI::SternheimerFDHamiltonian::Grid grid;
    grid.nx = 2;
    grid.ny = 1;
    grid.nz = 1;
    grid.hx = 1.0;
    grid.hy = 1.0;
    grid.hz = 1.0;
    grid.periodic = false;

    const std::vector<std::vector<ModuleRI::SternheimerRadialPerturbation>> radials_by_type = {{radial}};
    const std::vector<int> atom_types = {0};
    const std::vector<ModuleBase::Vector3<double>> atom_positions = {ModuleBase::Vector3<double>(0.0, 0.0, 0.0)};

    const std::vector<ModuleRI::SternheimerABFGridChannel> channels
        = ModuleRI::sample_sternheimer_abf_grid_channels(radials_by_type, atom_types, atom_positions, grid, 1);

    ASSERT_EQ(channels.size(), 1);
    EXPECT_EQ(channels[0].atom_index, 0);
    EXPECT_EQ(channels[0].atom_local_index, 0);
    EXPECT_EQ(channels[0].angular_momentum, 0);
    EXPECT_EQ(channels[0].magnetic_index, 0);
    EXPECT_EQ(channels[0].label, "toy_s");
    ASSERT_EQ(channels[0].potential_r.size(), 2);

    constexpr double y00 = 0.28209479177387814347;
    EXPECT_NEAR(channels[0].potential_r[0], 2.0 * y00, 1.0e-14);
    EXPECT_NEAR(channels[0].potential_r[1], 4.0 * y00, 1.0e-14);
    EXPECT_NEAR(channels[0].max_abs, 4.0 * y00, 1.0e-14);
}

TEST(SternheimerABFSPerturbation, ReadsABFSDiagnosticOnlyEnvironmentFlag)
{
    unsetenv("ABACUS_STERNHEIMER_FD_ST_ABFS_DIAG_ONLY");
    EXPECT_FALSE(ModuleRI::sternheimer_abfs_diag_only_enabled());

    setenv("ABACUS_STERNHEIMER_FD_ST_ABFS_DIAG_ONLY", "1", 1);
    EXPECT_TRUE(ModuleRI::sternheimer_abfs_diag_only_enabled());

    setenv("ABACUS_STERNHEIMER_FD_ST_ABFS_DIAG_ONLY", "false", 1);
    EXPECT_FALSE(ModuleRI::sternheimer_abfs_diag_only_enabled());

    unsetenv("ABACUS_STERNHEIMER_FD_ST_ABFS_DIAG_ONLY");
}
