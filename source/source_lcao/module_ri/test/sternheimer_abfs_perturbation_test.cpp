#include "source_lcao/module_ri/sternheimer_abfs_perturbation.h"
#include "source_lcao/module_ri/sternheimer_abacus_st_smoke.h"

#include <gtest/gtest.h>

#include <cmath>
#include <complex>
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

TEST(SternheimerABFSPerturbation, RequiresExplicitOptInForLargeGridCoulombDiagnostic)
{
    unsetenv("ABACUS_STERNHEIMER_GRID_COULOMB_DIAG");
    EXPECT_TRUE(ModuleRI::sternheimer_grid_coulomb_diagnostic_enabled(32));
    EXPECT_FALSE(ModuleRI::sternheimer_grid_coulomb_diagnostic_enabled(33));

    setenv("ABACUS_STERNHEIMER_GRID_COULOMB_DIAG", "1", 1);
    EXPECT_TRUE(ModuleRI::sternheimer_grid_coulomb_diagnostic_enabled(244));

    setenv("ABACUS_STERNHEIMER_GRID_COULOMB_DIAG", "false", 1);
    EXPECT_FALSE(ModuleRI::sternheimer_grid_coulomb_diagnostic_enabled(244));
    unsetenv("ABACUS_STERNHEIMER_GRID_COULOMB_DIAG");
}

TEST(SternheimerABFSPerturbation, BlochPotentialUsesQPhaseAndGammaMatchesRealSampler)
{
    ModuleRI::SternheimerRadialPerturbation radial;
    radial.angular_momentum = 0;
    radial.radial_index = 0;
    radial.label = "toy_s";
    radial.radial_grid = {0.0, 0.5, 1.0};
    radial.radial_values = {0.0, 1.0, 2.0};

    ModuleRI::SternheimerFDHamiltonian::Grid grid;
    grid.nx = 3;
    grid.ny = 1;
    grid.nz = 1;
    grid.hx = 1.0;
    grid.hy = 10.0;
    grid.hz = 10.0;
    grid.periodic = true;

    const std::vector<std::vector<ModuleRI::SternheimerRadialPerturbation>> radials_by_type = {{radial}};
    const std::vector<int> atom_types = {0};
    const std::vector<ModuleBase::Vector3<double>> atom_positions
        = {ModuleBase::Vector3<double>(0.0, 0.0, 0.0)};
    const auto bloch = ModuleRI::sample_sternheimer_abf_bloch_grid_channels(
        radials_by_type, atom_types, atom_positions, grid, {0.25, 0.0, 0.0}, 1);

    constexpr double y00 = 0.28209479177387814347;
    ASSERT_EQ(bloch.size(), 1);
    ASSERT_EQ(bloch[0].potential_r.size(), 3);
    EXPECT_NEAR(bloch[0].potential_r[2].real(), 0.0, 1.0e-14);
    EXPECT_NEAR(bloch[0].potential_r[2].imag(), 2.0 * y00, 1.0e-14);

    const auto gamma_complex = ModuleRI::sample_sternheimer_abf_bloch_grid_channels(
        radials_by_type, atom_types, atom_positions, grid, {0.0, 0.0, 0.0}, 1);
    const auto gamma_real
        = ModuleRI::sample_sternheimer_abf_grid_channels(radials_by_type, atom_types, atom_positions, grid, 1);
    ASSERT_EQ(gamma_complex.size(), gamma_real.size());
    for (std::size_t ir = 0; ir != gamma_real[0].potential_r.size(); ++ir)
    {
        EXPECT_NEAR(gamma_complex[0].potential_r[ir].real(), gamma_real[0].potential_r[ir], 1.0e-14);
        EXPECT_NEAR(gamma_complex[0].potential_r[ir].imag(), 0.0, 1.0e-14);
    }
}

TEST(SternheimerABFSPerturbation, PeriodicPoissonAppliesShiftedFullCoulombKernel)
{
    ModuleRI::SternheimerFDHamiltonian::Grid grid;
    grid.nx = 4;
    grid.ny = 2;
    grid.nz = 2;
    grid.hx = 0.5 * M_PI;
    grid.hy = M_PI;
    grid.hz = M_PI;
    grid.periodic = true;

    ModuleRI::SternheimerABFBlochGridChannel density;
    density.channel_index = 0;
    density.atom_index = 0;
    density.atom_local_index = 0;
    density.angular_momentum = 0;
    density.label = "plane_wave";
    const std::size_t grid_size
        = static_cast<std::size_t>(grid.nx) * static_cast<std::size_t>(grid.ny)
          * static_cast<std::size_t>(grid.nz);
    density.potential_r.resize(grid_size);

    const ModuleRI::SternheimerReducedKPoint qpoint{0.5, 0.0, 0.0};
    for (int ix = 0; ix != grid.nx; ++ix)
    {
        const double x = ix * grid.hx;
        for (int iy = 0; iy != grid.ny; ++iy)
        {
            for (int iz = 0; iz != grid.nz; ++iz)
            {
                const std::size_t ir
                    = static_cast<std::size_t>((ix * grid.ny + iy) * grid.nz + iz);
                density.potential_r[ir] = std::exp(std::complex<double>(0.0, 1.5 * x));
            }
        }
    }

    const auto potentials
        = ModuleRI::solve_sternheimer_abf_periodic_full_coulomb({density}, grid, qpoint);
    ASSERT_EQ(potentials.size(), 1);
    ASSERT_EQ(potentials[0].potential_r.size(), density.potential_r.size());
    const double expected_factor = 4.0 * M_PI / (1.5 * 1.5);
    for (std::size_t ir = 0; ir != density.potential_r.size(); ++ir)
    {
        EXPECT_NEAR(potentials[0].potential_r[ir].real(),
                    expected_factor * density.potential_r[ir].real(),
                    1.0e-12);
        EXPECT_NEAR(potentials[0].potential_r[ir].imag(),
                    expected_factor * density.potential_r[ir].imag(),
                    1.0e-12);
    }
}

TEST(SternheimerABFSPerturbation, NonOrthogonalSamplingUsesCartesianLatticeCombination)
{
    ModuleRI::SternheimerRadialPerturbation radial;
    radial.angular_momentum = 0;
    radial.radial_index = 0;
    radial.label = "linear_s";
    radial.radial_grid = {0.0, 1.0};
    radial.radial_values = {0.0, 1.0};

    ModuleRI::SternheimerFDHamiltonian::Grid grid;
    grid.nx = 4;
    grid.ny = 4;
    grid.nz = 4;
    grid.hx = std::sqrt(2.0) / 4.0;
    grid.hy = std::sqrt(2.0) / 4.0;
    grid.hz = std::sqrt(2.0) / 4.0;
    grid.periodic = false;
    grid.lattice_vectors = {{{0.0, 1.0, 1.0},
                             {1.0, 0.0, 1.0},
                             {1.0, 1.0, 0.0}}};

    const std::vector<std::vector<ModuleRI::SternheimerRadialPerturbation>> radials_by_type = {{radial}};
    const std::vector<int> atom_types = {0};
    const std::vector<ModuleBase::Vector3<double>> atom_positions
        = {ModuleBase::Vector3<double>(0.0, 0.0, 0.0)};
    const auto channels
        = ModuleRI::sample_sternheimer_abf_grid_channels(radials_by_type, atom_types, atom_positions, grid, 1);

    ASSERT_EQ(channels.size(), 1);
    const std::size_t ir = static_cast<std::size_t>((1 * grid.ny + 1) * grid.nz);
    constexpr double y00 = 0.28209479177387814347;
    EXPECT_NEAR(channels[0].potential_r[ir], std::sqrt(0.375) * y00, 1.0e-14);
}

TEST(SternheimerABFSPerturbation, NonOrthogonalPeriodicPoissonUsesReciprocalMetric)
{
    ModuleRI::SternheimerFDHamiltonian::Grid grid;
    grid.nx = 4;
    grid.ny = 4;
    grid.nz = 4;
    grid.hx = std::sqrt(2.0) / 4.0;
    grid.hy = std::sqrt(2.0) / 4.0;
    grid.hz = std::sqrt(2.0) / 4.0;
    grid.periodic = true;
    grid.lattice_vectors = {{{0.0, 1.0, 1.0},
                             {1.0, 0.0, 1.0},
                             {1.0, 1.0, 0.0}}};

    ModuleRI::SternheimerABFBlochGridChannel density;
    density.channel_index = 0;
    density.atom_index = 0;
    density.atom_local_index = 0;
    density.angular_momentum = 0;
    density.label = "skew_plane_wave";
    density.potential_r.resize(static_cast<std::size_t>(grid.nx * grid.ny * grid.nz));

    const ModuleRI::SternheimerReducedKPoint qpoint{0.25, 0.0, 0.0};
    for (int ix = 0; ix != grid.nx; ++ix)
    {
        for (int iy = 0; iy != grid.ny; ++iy)
        {
            for (int iz = 0; iz != grid.nz; ++iz)
            {
                const double reduced_phase
                    = 1.25 * static_cast<double>(ix) / grid.nx
                      + static_cast<double>(iy) / grid.ny;
                const std::size_t ir = static_cast<std::size_t>((ix * grid.ny + iy) * grid.nz + iz);
                density.potential_r[ir]
                    = std::exp(std::complex<double>(0.0, 2.0 * M_PI * reduced_phase));
            }
        }
    }

    const auto potentials
        = ModuleRI::solve_sternheimer_abf_periodic_full_coulomb({density}, grid, qpoint);
    ASSERT_EQ(potentials.size(), 1);
    const double wavevector_squared = 4.0 * M_PI * M_PI * 1.296875;
    const double expected_factor = 4.0 * M_PI / wavevector_squared;
    for (std::size_t ir = 0; ir != density.potential_r.size(); ++ir)
    {
        EXPECT_NEAR(potentials[0].potential_r[ir].real(),
                    expected_factor * density.potential_r[ir].real(),
                    1.0e-12);
        EXPECT_NEAR(potentials[0].potential_r[ir].imag(),
                    expected_factor * density.potential_r[ir].imag(),
                    1.0e-12);
    }
}

TEST(SternheimerABFSPerturbation, PeriodicPoissonRejectsGammaBodyKernel)
{
    ModuleRI::SternheimerFDHamiltonian::Grid grid;
    grid.nx = 2;
    grid.ny = 2;
    grid.nz = 2;
    grid.hx = 1.0;
    grid.hy = 1.0;
    grid.hz = 1.0;
    grid.periodic = true;

    ModuleRI::SternheimerABFBlochGridChannel density;
    const std::size_t grid_size
        = static_cast<std::size_t>(grid.nx) * static_cast<std::size_t>(grid.ny)
          * static_cast<std::size_t>(grid.nz);
    density.potential_r.assign(grid_size, {1.0, 0.0});
    EXPECT_THROW(ModuleRI::solve_sternheimer_abf_periodic_full_coulomb(
                     {density}, grid, {0.0, 0.0, 0.0}),
                 std::invalid_argument);
}
