#include "../kedf_vw.h"

#include "source_base/constants.h"

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

class TestKEDFVW : public ::testing::Test
{
protected:
    void SetUp() override
    {
        // A cubic cell of side 2*pi makes the fundamental z wave number one.
        const double length = ModuleBase::TWO_PI;
        const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
#ifdef __MPI
        basis.initmpi(1, 0, MPI_COMM_SELF);
#endif
        basis.initgrids(length, lattice, 8, 8, 8);
        basis.initparameters(false, 4.0, 1, false);
        basis.setuptransform();
        basis.collect_local_pw();
        volume = length * length * length;
        const double dV = volume / basis.nxyz;
        kedf.set_para(dV, weight);

        phi.resize(basis.nrxx);
        cosine.resize(basis.nrxx);
        for (int ir = 0; ir < basis.nrxx; ++ir)
        {
            const int iz = ir % basis.nplane + basis.startz_current;
            const double phase = ModuleBase::TWO_PI * iz / basis.nz;
            cosine[ir] = std::cos(phase);
            phi[ir] = offset + amplitude * cosine[ir];
        }
    }

    ModulePW::PW_Basis basis{"cpu", "double"};
    KEDF_vW kedf;
    std::vector<double> phi;
    std::vector<double> cosine;
    double volume = 0.0;
    const double offset = 2.0;
    const double amplitude = 0.5;
    const double weight = 0.75;
};

TEST_F(TestKEDFVW, PeriodicEnergyMatchesAnalyticIntegral)
{
    // In Ry, E = weight * integral(phi * (-laplacian(phi)))
    //          = weight * volume * amplitude^2 / 2 for wave number one.
    double* input[] = {phi.data()};
    const double expected = weight * volume * amplitude * amplitude / 2.0;
    EXPECT_NEAR(kedf.get_energy(input, &basis), expected, 1e-10);
    EXPECT_NEAR(kedf.vw_energy, expected, 1e-10);

    for (double& value : phi)
    {
        value = -value;
    }
    EXPECT_NEAR(kedf.get_energy(input, &basis), expected, 1e-10);
}

TEST_F(TestKEDFVW, PeriodicEnergyDensityMatchesNegativeLaplacian)
{
    double* input[] = {phi.data()};
    const int ir = 0;
    const double expected = weight * phi[ir] * amplitude * cosine[ir];
    EXPECT_NEAR(kedf.get_energy_density(input, 0, ir, &basis), expected, 1e-10);

    for (double& value : phi)
    {
        value = -value;
    }
    EXPECT_NEAR(kedf.get_energy_density(input, 0, ir, &basis), expected, 1e-10);
}

TEST_F(TestKEDFVW, PeriodicPotentialAccumulatesAndChangesSign)
{
    const double* input[] = {phi.data()};
    ModuleBase::matrix potential(1, basis.nrxx);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        potential(0, ir) = 1.25;
    }
    kedf.vw_potential(input, &basis, potential);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const double expected = 1.25 + 2.0 * weight * amplitude * cosine[ir];
        EXPECT_NEAR(potential(0, ir), expected, 1e-10);
        phi[ir] = -phi[ir];
    }

    // The negative input contributes the opposite potential, cancelling it.
    kedf.vw_potential(input, &basis, potential);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        EXPECT_NEAR(potential(0, ir), 1.25, 1e-10);
    }
    const double expected_energy = weight * volume * amplitude * amplitude / 2.0;
    EXPECT_NEAR(kedf.vw_energy, expected_energy, 1e-10);
}

TEST_F(TestKEDFVW, PeriodicStressHasOnlyLongitudinalComponent)
{
    const double* input[] = {phi.data()};
    const double expected = weight * amplitude * amplitude;
    kedf.get_stress(input, &basis);
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            if (i == 2 && j == 2)
            {
                EXPECT_NEAR(kedf.stress(i, j), expected, 1e-10);
            }
            else
            {
                EXPECT_NEAR(kedf.stress(i, j), 0.0, 1e-10);
            }
        }
    }
    for (double& value : phi)
    {
        value = -value;
    }
    kedf.get_stress(input, &basis);
    EXPECT_NEAR(kedf.stress(2, 2), expected, 1e-10);
}
