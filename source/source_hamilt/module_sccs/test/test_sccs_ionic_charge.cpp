#include "sccs_test.h"
#include "../sccs_ionic_charge.h"

#include "source_base/parallel_reduce.h"
#include "source_cell/cell_tools.h"

using SccsIonicChargeTest = SccsTest::PwTest;

TEST_F(SccsIonicChargeTest, ValenceNormalizationAndFourierPhase)
{
    std::vector<unitcell::AtomData> atoms(2);
    atoms[0].position.x = length / 4.0;
    atoms[0].valence_charge = 8.0;
    atoms[1].valence_charge = 1.0;
    const double spread = ModuleSccs::gaussian_ion_spread;
    std::vector<double> density;
    ModuleSccs::gaussian_ionic_density(atoms, basis, tpiba, spread, density);
    double integral = 0.0;
    for (double value : density)
    {
        integral += value;
    }
    Parallel_Reduce::reduce_pool(integral);
    integral *= basis.omega / basis.nxyz;
    EXPECT_NEAR(integral, 9.0, 1e-12);

    std::vector<std::complex<double>> charge_g(basis.npw);
    basis.real2recip(density.data(), charge_g.data());
    const double exponent = -spread * spread * tpiba * tpiba / 4.0;
    const double amplitude = 8.0 / basis.omega * std::exp(exponent);
    const double second_amplitude = 1.0 / basis.omega * std::exp(exponent);
    double selected = 0.0;
    for (int ig = 0; ig < basis.npw; ++ig)
    {
        const ModuleBase::Vector3<double>& g = basis.gdirect[ig];
        if (g.x == 1.0 && g.y == 0.0 && g.z == 0.0)
        {
            selected = 1.0;
            EXPECT_NEAR(charge_g[ig].real(), second_amplitude, 1e-14);
            const double expected_imaginary = -amplitude;
            EXPECT_NEAR(charge_g[ig].imag(), expected_imaginary, 1e-14);
        }
    }
    Parallel_Reduce::reduce_pool(selected);
    EXPECT_DOUBLE_EQ(selected, 1.0);
    atoms[0].position.x += length;
    std::vector<double> translated;
    ModuleSccs::gaussian_ionic_density(atoms, basis, tpiba, spread, translated);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        EXPECT_NEAR(translated[ir], density[ir], 1e-13);
    }
}
