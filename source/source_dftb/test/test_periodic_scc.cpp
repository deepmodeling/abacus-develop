#include "source_dftb/periodic_scc.h"

#include "gtest/gtest.h"

#include <array>
#include <vector>

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

} // namespace ModuleDFTB
