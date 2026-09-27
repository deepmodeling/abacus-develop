#include "source_dftb/eigensolver.h"

#include "gtest/gtest.h"

#include <cmath>
#include <complex>
#include <vector>

namespace ModuleDFTB
{

TEST(DftbNativeEigensolverTest, SolvesGeneralizedHermitianProblem)
{
    DftbBlochMatrices matrices;
    matrices.dimension = 2;
    matrices.hamiltonian.resize(4, std::complex<double>(0.0, 0.0));
    matrices.overlap.resize(4, std::complex<double>(0.0, 0.0));
    matrices.h(0, 0) = 1.0;
    matrices.h(1, 1) = 4.0;
    matrices.s(0, 0) = 1.0;
    matrices.s(1, 1) = 2.0;

    const DftbEigenSolution solution = solve_generalized_hermitian(matrices);

    ASSERT_EQ(solution.eigenvalues_hartree.size(), 2U);
    EXPECT_NEAR(solution.eigenvalues_hartree[0], 1.0, 1.0e-12);
    EXPECT_NEAR(solution.eigenvalues_hartree[1], 2.0, 1.0e-12);
}

TEST(DftbNativeEigensolverTest, FiniteTemperatureFillingConservesElectrons)
{
    std::vector<DftbKPointSpectrum> spectra(2);
    spectra[0].weight = 0.5;
    spectra[0].solution.dimension = 2;
    spectra[0].solution.eigenvalues_hartree = {0.0, 1.0};
    spectra[1].weight = 0.5;
    spectra[1].solution.dimension = 2;
    spectra[1].solution.eigenvalues_hartree = {0.2, 1.2};

    const DftbFermiFilling filling = fill_fermi_occupations(spectra, 2.0, 0.01);

    double weighted_electrons = 0.0;
    for (std::size_t k = 0; k < spectra.size(); ++k)
    {
        for (std::size_t band = 0; band < spectra[k].solution.dimension; ++band)
        {
            const double occupation = filling.occupations[k][band];
            EXPECT_GE(occupation, 0.0);
            EXPECT_LE(occupation, 2.0);
            weighted_electrons += spectra[k].weight * occupation;
        }
    }
    EXPECT_NEAR(weighted_electrons, 2.0, 1.0e-10);
}

TEST(DftbNativeEigensolverTest, ZeroWeightKPointsDoNotAffectFilling)
{
    std::vector<DftbKPointSpectrum> spectra(2);
    spectra[0].weight = 0.0;
    spectra[0].solution.dimension = 2;
    spectra[0].solution.eigenvalues_hartree = {-10.0, -9.0};
    spectra[1].weight = 1.0;
    spectra[1].solution.dimension = 2;
    spectra[1].solution.eigenvalues_hartree = {0.0, 1.0};

    const DftbFermiFilling filling = fill_fermi_occupations(spectra, 2.0, 0.0);

    EXPECT_NEAR(filling.fermi_energy_hartree, 0.0, 1.0e-12);
    EXPECT_DOUBLE_EQ(filling.occupations[0][0], 0.0);
    EXPECT_DOUBLE_EQ(filling.occupations[0][1], 0.0);
    EXPECT_DOUBLE_EQ(filling.occupations[1][0], 2.0);
    EXPECT_DOUBLE_EQ(filling.occupations[1][1], 0.0);
    EXPECT_NEAR(filling.band_energy_hartree, 0.0, 1.0e-12);
}

} // namespace ModuleDFTB
