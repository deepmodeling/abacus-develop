#include "../ovlp_block.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_basis/module_nao/two_center_integrator.h"
#include "source_cell/unitcell.h"
#include <gtest/gtest.h>
#include <complex>
#include <memory>

class OverlapBlockTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        atom.reset(new Atom);
        cell.ntype = 1;
        cell.nat = 1;
        cell.atoms = atom.get();
        // UnitCell owns the atom-index maps as vectors.
        cell.iat2it = {0};
        cell.iat2ia = {0};
        atom->nw = 2;
        atom->iw2l = {0, 0};
        atom->iw2n = {0, 0};
        atom->iw2m = {0, 0};
    }
    void set_polarization(int npol)
    {
        const std::vector<int> offsets{0};
        cell.set_iat2iwt_for_test(offsets, npol);
        const int size = atom->nw * npol;
        distribution.set_serial(size, size);
        distribution.set_atomic_trace(cell.get_iat2iwt(), cell.nat, size);
    }
    UnitCell cell;
    std::unique_ptr<Atom> atom;
    Parallel_Orbitals distribution;
    TwoCenterIntegrator integrator;
};

TEST_F(OverlapBlockTest, AddsToExistingRealBlock)
{
    set_polarization(1);
    const ModuleBase::Vector3<double> displacement(0.3, -0.2, 0.1);
    std::vector<double> values(4, 2.0);
    // The existing operator fixture mocks each two-center integral as one.
    hamilt::cal_overlap_block(cell, integrator, 0, 0, distribution, displacement, values.data());
    hamilt::cal_overlap_block(cell, integrator, 0, 0, distribution, displacement, values.data());
    for (double value : values) { EXPECT_DOUBLE_EQ(value, 4.0); }
}

TEST_F(OverlapBlockTest, PreservesComplexSpinOffDiagonalEntries)
{
    set_polarization(2);
    const ModuleBase::Vector3<double> displacement(-0.3, 0.2, -0.1);
    const std::complex<double> initial(2.0, 3.0);
    std::vector<std::complex<double>> values(16, initial);
    hamilt::cal_overlap_block(cell, integrator, 0, 0, distribution, displacement, values.data());
    for (int row = 0; row < 4; ++row)
    {
        for (int col = 0; col < 4; ++col)
        {
            const double increment = row % 2 == col % 2 ? 1.0 : 0.0;
            const std::complex<double> expected(2.0 + increment, 3.0);
            EXPECT_EQ(values[row * 4 + col], expected);
        }
    }
}
