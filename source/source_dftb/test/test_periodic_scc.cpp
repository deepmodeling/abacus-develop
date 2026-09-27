#include "source_dftb/periodic_scc.h"

#include "gtest/gtest.h"

#include <array>
#include <limits>
#include <stdexcept>
#include <vector>

#ifdef __MPI
#include <mpi.h>

namespace
{
class DftbMpiTestEnvironment : public ::testing::Environment
{
  public:
    void SetUp() override
    {
        int initialized = 0;
        MPI_Initialized(&initialized);
        if (!initialized)
        {
            int argc = 0;
            char** argv = nullptr;
            MPI_Init(&argc, &argv);
            this->owns_mpi_ = true;
        }
    }

    void TearDown() override
    {
        int finalized = 0;
        MPI_Finalized(&finalized);
        if (this->owns_mpi_ && !finalized) MPI_Finalize();
    }

  private:
    bool owns_mpi_ = false;
};

::testing::Environment* const dftb_mpi_test_environment =
    ::testing::AddGlobalTestEnvironment(new DftbMpiTestEnvironment);
} // namespace
#endif

namespace ModuleDFTB
{

TEST(DftbNativePeriodicSccTest, SolvesNeutralCellAndFrozenPotentialBands)
{
    SkfData hydrogen_like;
    hydrogen_like.filename = "synthetic-H.skf";
    hydrogen_like.homonuclear = true;
    hydrogen_like.grid_spacing_bohr = 0.1;
    hydrogen_like.declared_grid_points = 1;
    hydrogen_like.has_atomic_data = true;
    hydrogen_like.onsite_hartree = {{-0.5, 0.5, 0.0}};
    hydrogen_like.hubbard_u_hartree = {{0.4, 0.4, 0.4}};
    // The reader stores SKF homonuclear occupations in s, p, d order.
    hydrogen_like.reference_occupations = {{2.0, 0.0, 0.0}};

    DftbPeriodicInput input;
    input.atoms.resize(1);
    input.atoms[0].species = 0;
    input.atoms[0].homonuclear_data = &hydrogen_like;
    DftbPairParameters pair;
    pair.species_a = 0;
    pair.species_b = 0;
    pair.ab = &hydrogen_like;
    pair.ba = &hydrogen_like;
    input.pair_parameters.push_back(pair);
    input.lattice_bohr[0] = {{40.0, 0.0, 0.0}};
    input.lattice_bohr[1] = {{0.0, 40.0, 0.0}};
    input.lattice_bohr[2] = {{0.0, 0.0, 40.0}};
    DftbWeightedKPoint gamma;
    gamma.fractional = {{0.0, 0.0, 0.0}};
    gamma.weight = 1.0;
    input.kpoints.push_back(gamma);
    DftbBandKPoint band_gamma;
    band_gamma.fractional = {{0.0, 0.0, 0.0}};
    band_gamma.label = "Gamma";
    input.band_kpoints.push_back(band_gamma);
    DftbBandKPoint band_x;
    band_x.fractional = {{0.5, 0.0, 0.0}};
    band_x.label = "X";
    input.band_kpoints.push_back(band_x);
    input.hubbard_derivative.push_back(-0.1);
    input.total_electrons = 2.0;
    input.maximum_scc_iterations = 10;
    input.scc_tolerance = 1.0e-10;
    input.third_order = true;

    const DftbPeriodicResult result = solve_periodic_dftb(input);

    ASSERT_TRUE(result.converged);
    ASSERT_EQ(result.scc_iterations, 1);
    ASSERT_EQ(result.kpoint_eigenvalues.size(), 1U);
    ASSERT_EQ(result.kpoint_eigenvalues[0].eigenvalues_hartree.size(), 4U);
    EXPECT_NEAR(result.electron_excess_charges[0], 0.0, 1.0e-12);
    EXPECT_NEAR(result.total_free_energy_hartree, -1.0, 1.0e-10);
    ASSERT_EQ(result.band_structure.size(), 2U);
    ASSERT_EQ(result.band_structure[0].eigenvalues_hartree.size(), 4U);
    EXPECT_NEAR(result.band_structure[0].eigenvalues_hartree[0], -0.5, 1.0e-12);
    EXPECT_NEAR(result.band_structure[1].eigenvalues_hartree[0], -0.5, 1.0e-12);
    EXPECT_NEAR(result.band_structure[1].distance_inverse_bohr, 0.07853981633974483, 1.0e-12);
}

TEST(DftbNativePeriodicSccTest, IterationEnergyUsesThePotentialThatGeneratedItsDensity)
{
    const auto initialize_skf = [](SkfData* data, const std::string& filename, bool homonuclear) {
        data->filename = filename;
        data->homonuclear = homonuclear;
        data->grid_spacing_bohr = 0.1;
        data->declared_grid_points = 9;
        data->hamiltonian.resize(8);
        data->overlap.resize(8);
        data->has_repulsive_spline = true;
        data->repulsive.cutoff_bohr = 1.0;
        data->repulsive.interval_starts_bohr = {0.0, 0.5};
        data->repulsive.interval_ends_bohr = {0.5, 1.0};
        data->repulsive.cubic_coefficients.resize(1);
        if (homonuclear)
        {
            data->has_atomic_data = true;
            data->hubbard_u_hartree = {{0.4, 0.4, 0.4}};
            data->reference_occupations = {{1.0, 0.0, 0.0}};
        }
    };

    SkfData species_a;
    SkfData species_b;
    SkfData a_to_b;
    SkfData b_to_a;
    initialize_skf(&species_a, "synthetic-A-A.skf", true);
    initialize_skf(&species_b, "synthetic-B-B.skf", true);
    initialize_skf(&a_to_b, "synthetic-A-B.skf", false);
    initialize_skf(&b_to_a, "synthetic-B-A.skf", false);
    species_a.onsite_hartree = {{-1000.0, -999.0, 0.0}};
    species_b.onsite_hartree = {{1000.0, 1001.0, 0.0}};

    DftbPeriodicInput input;
    input.atoms.resize(2);
    input.atoms[0].species = 0;
    input.atoms[0].homonuclear_data = &species_a;
    input.atoms[1].species = 1;
    input.atoms[1].homonuclear_data = &species_b;
    input.atoms[1].position_bohr = {{2.0, 0.0, 0.0}};
    const auto add_pair = [&input](std::size_t a, std::size_t b, const SkfData* ab, const SkfData* ba) {
        DftbPairParameters pair;
        pair.species_a = a;
        pair.species_b = b;
        pair.ab = ab;
        pair.ba = ba;
        input.pair_parameters.push_back(pair);
    };
    add_pair(0, 0, &species_a, &species_a);
    add_pair(0, 1, &a_to_b, &b_to_a);
    add_pair(1, 1, &species_b, &species_b);
    input.lattice_bohr[0] = {{40.0, 0.0, 0.0}};
    input.lattice_bohr[1] = {{0.0, 40.0, 0.0}};
    input.lattice_bohr[2] = {{0.0, 0.0, 40.0}};
    DftbWeightedKPoint gamma;
    gamma.weight = 1.0;
    input.kpoints.push_back(gamma);
    input.hubbard_derivative = {0.0, 0.0};
    input.total_electrons = 2.0;
    input.maximum_scc_iterations = 4;
    input.scc_tolerance = 1.0e-10;
    input.mixing_parameter = 1.0;

    DftbPeriodicInput invalid_input = input;
    invalid_input.scc_tolerance = std::numeric_limits<double>::infinity();
    EXPECT_THROW(solve_periodic_dftb(invalid_input), std::invalid_argument);

    std::vector<DftbSccIteration> iterations;
    const DftbPeriodicResult result = solve_periodic_dftb(
        input, [&iterations](const DftbSccIteration& iteration) { iterations.push_back(iteration); });

    ASSERT_TRUE(result.converged);
    ASSERT_EQ(iterations.size(), 2U);
    EXPECT_GT(iterations.front().maximum_charge_residual, input.scc_tolerance);
    EXPECT_NEAR(iterations.front().electronic_energy_hartree,
                result.h0_energy_hartree + result.scc_energy_hartree,
                1.0e-10);
}

TEST(DftbNativePeriodicSccTest, RejectsOverlappingDistinctAtoms)
{
    SkfData hydrogen_like;
    hydrogen_like.filename = "synthetic-H.skf";
    hydrogen_like.homonuclear = true;
    hydrogen_like.grid_spacing_bohr = 0.1;
    hydrogen_like.declared_grid_points = 1;
    hydrogen_like.has_atomic_data = true;
    hydrogen_like.onsite_hartree = {{-0.5, 0.5, 0.0}};
    hydrogen_like.hubbard_u_hartree = {{0.4, 0.4, 0.4}};
    hydrogen_like.reference_occupations = {{2.0, 0.0, 0.0}};

    DftbPeriodicInput input;
    input.atoms.resize(2);
    for (auto& atom : input.atoms)
    {
        atom.species = 0;
        atom.homonuclear_data = &hydrogen_like;
    }
    DftbPairParameters pair;
    pair.species_a = 0;
    pair.species_b = 0;
    pair.ab = &hydrogen_like;
    pair.ba = &hydrogen_like;
    input.pair_parameters.push_back(pair);
    input.lattice_bohr[0] = {{40.0, 0.0, 0.0}};
    input.lattice_bohr[1] = {{0.0, 40.0, 0.0}};
    input.lattice_bohr[2] = {{0.0, 0.0, 40.0}};
    DftbWeightedKPoint gamma;
    gamma.fractional = {{0.0, 0.0, 0.0}};
    gamma.weight = 1.0;
    input.kpoints.push_back(gamma);
    input.hubbard_derivative.push_back(-0.1);
    input.total_electrons = 4.0;

    EXPECT_THROW(solve_periodic_dftb(input), std::runtime_error);
}

TEST(DftbNativePeriodicSccTest, RejectsUnsupportedDValenceOccupation)
{
    SkfData d_shell_species;
    d_shell_species.filename = "synthetic-spd.skf";
    d_shell_species.homonuclear = true;
    d_shell_species.has_atomic_data = true;
    d_shell_species.reference_occupations = {{1.0, 0.0, 1.0}};

    DftbPeriodicInput input;
    input.atoms.resize(1);
    input.atoms[0].species = 0;
    input.atoms[0].homonuclear_data = &d_shell_species;
    DftbPairParameters pair;
    pair.species_a = 0;
    pair.species_b = 0;
    pair.ab = &d_shell_species;
    pair.ba = &d_shell_species;
    input.pair_parameters.push_back(pair);
    input.hubbard_derivative.push_back(0.0);
    DftbWeightedKPoint gamma;
    gamma.fractional = {{0.0, 0.0, 0.0}};
    gamma.weight = 1.0;
    input.kpoints.push_back(gamma);
    input.total_electrons = 2.0;

    EXPECT_THROW(solve_periodic_dftb(input), std::invalid_argument);
}

} // namespace ModuleDFTB
