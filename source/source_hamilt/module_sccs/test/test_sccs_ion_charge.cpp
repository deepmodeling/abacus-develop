#include "sccs_test.h"
#include "../sccs_ion_charge.h"

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
    ModuleSccs::gauss_ion_rho(atoms, basis, tpiba, spread, density);
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
    ModuleSccs::gauss_ion_rho(atoms, basis, tpiba, spread, translated);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        EXPECT_NEAR(translated[ir], density[ir], 1e-13);
    }
}

// Atoms without a known atomic number keep the singleton width; a known Z ==
// zv atom is skipped by a singleton but not by a per-atom list; widths <= 0
// disable an atom; the atom data are never modified.
TEST_F(SccsIonicChargeTest, CoreGaussianSelection)
{
    std::vector<unitcell::AtomData> atoms(3);
    atoms[0].valence_charge = 6.0;
    atoms[1].valence_charge = 1.0;
    atoms[2].valence_charge = 7.0;
    atoms[0].position.x = length / 4.0;
    std::vector<double> widths = {0.5};
    std::vector<double> density;
    auto integral = [&]() {
        double sum = 0.0;
        for (double value : density)
        {
            sum += value;
        }
        Parallel_Reduce::reduce_pool(sum);
        return sum * basis.omega / basis.nxyz;
    };
    ModuleSccs::gauss_core_rho(atoms, basis, tpiba, widths, density);
    EXPECT_NEAR(integral(), 14.0, 1e-12);
    widths = {0.8};
    ModuleSccs::gauss_core_rho(atoms, basis, tpiba, widths, density);
    EXPECT_NEAR(integral(), 14.0, 1e-12);
    widths = {0.7, 0.0, -1.0};
    ModuleSccs::gauss_core_rho(atoms, basis, tpiba, widths, density);
    EXPECT_NEAR(integral(), 6.0, 1e-12);
    std::vector<unitcell::AtomData> selected;
    std::vector<double> resolved;
    ModuleSccs::prep_core_gauss(atoms, widths, selected, resolved);
    EXPECT_DOUBLE_EQ(resolved[0], 0.7);
    EXPECT_DOUBLE_EQ(selected[1].valence_charge, 0.0);
    EXPECT_DOUBLE_EQ(atoms[1].valence_charge, 1.0);
    widths = {0.0};
    ModuleSccs::gauss_core_rho(atoms, basis, tpiba, widths, density);
    EXPECT_DOUBLE_EQ(integral(), 0.0);

    std::vector<unitcell::AtomData> known(4);
    known[0].atomic_number = 8;
    known[0].valence_charge = 6.0;
    known[1].atomic_number = 1;
    known[1].valence_charge = 1.0;
    known[2].atomic_number = 2;
    known[2].valence_charge = 2.0;
    known[3].valence_charge = 3.0; // Unknown element: do not silently exclude.
    widths = {0.8};
    ModuleSccs::gauss_core_rho(known, basis, tpiba, widths, density);
    EXPECT_NEAR(integral(), 9.0, 1e-12);
    widths = {0.8, 0.8, 0.8, 0.8};
    ModuleSccs::gauss_core_rho(known, basis, tpiba, widths, density);
    EXPECT_NEAR(integral(), 12.0, 1e-12);
    widths = {0.8, 0.0, -1.0, 0.8};
    ModuleSccs::gauss_core_rho(known, basis, tpiba, widths, density);
    EXPECT_NEAR(integral(), 9.0, 1e-12);
    EXPECT_DOUBLE_EQ(known[1].valence_charge, 1.0);
}
