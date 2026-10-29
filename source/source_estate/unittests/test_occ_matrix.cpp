#include "source_estate/occ_matrix.h"

#include "source_cell/atom_spec.h"
#include "source_cell/magnetism.h"
#include "source_cell/unitcell.h"
#include "gtest/gtest.h"

#include <vector>

// UnitCell's constructor/destructor reference Magnetism; provide the
// minimal mock definitions (the real magnetism.cpp pulls heavy deps).
Magnetism::Magnetism()
{
    this->tot_mag = 0.0;
    this->abs_mag = 0.0;
}
Magnetism::~Magnetism()
{
}

/***********************************************************************
 * Unit tests for OccupationMatrix.
 *
 * Covered:
 * 1. init() layout for nspin=1, 2, 4 (channel count and matrix dims).
 * 2. init() skips uncorrelated atom types (l_channel=-1).
 * 3. init() builds corr_iwt with atom offsets and spinor strides from
 *    the first radial channel (n=0) only.
 * 4. write_save_to_flat() for nspin=1/2/4, including the empty-buffer
 *    no-op path.
 *
 * A minimal UnitCell is built directly in each fixture; no private
 * helpers or access-control tricks are used.
 ***********************************************************************/

class OccMatrixTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
    }

    /// Build a minimal UnitCell with `ntype` types. Each type has one
    /// atom. `l_channels` gives the correlated l for each type (-1 for
    /// uncorrelated). `nwl` and `l_nchi` are set to match `l_channels`.
    void make_cell(const std::vector<int>& l_channels,
                   const std::vector<int>& nwls,
                   const std::vector<std::vector<int>>& l_nchis,
                   int npol)
    {
        cell.ntype = static_cast<int>(l_channels.size());
        cell.nat = cell.ntype;
        atoms_storage.resize(cell.ntype);
        cell.atoms = atoms_storage.data();
        cell.iat2it.resize(cell.nat);
        cell.iat2ia.resize(cell.nat);
        cell.itia2iat.create(cell.ntype, cell.nat);

        for (int it = 0; it < cell.ntype; it++)
        {
            cell.atoms[it].na = 1;
            cell.atoms[it].nwl = nwls[it];
            cell.atoms[it].l_nchi = l_nchis[it];
            cell.atoms[it].nw = 0;
            for (int l = 0; l <= nwls[it]; l++)
            {
                cell.atoms[it].nw += (2 * l + 1) * l_nchis[it][l];
            }
            cell.atoms[it].set_index();
            cell.iat2it[it] = it;
            cell.iat2ia[it] = 0;
            cell.itia2iat(it, 0) = it;
        }
        cell.set_iat2iwt(npol);
    }

    UnitCell cell;
    std::vector<Atom> atoms_storage;
};

// ----------------------------------------------------------------------
// init: nspin=1 layout
// ----------------------------------------------------------------------
TEST_F(OccMatrixTest, InitCollinearNspin1)
{
    make_cell({2}, {2}, {{0, 0, 1}}, 1);
    OccupationMatrix occmat;
    occmat.init(cell, {2}, 1, 1);

    EXPECT_EQ(occmat.nspin(), 1);
    EXPECT_EQ(occmat.npol(), 1);

    const int iat = 0;
    const int l = 2;
    EXPECT_EQ(occmat.data()[iat][l].size(), 2u);
    EXPECT_EQ(occmat.data()[iat][l][0].nr, 5);
    EXPECT_EQ(occmat.data()[iat][l][0].nc, 5);
    EXPECT_EQ(occmat.data()[iat][l][1].nr, 5);
    EXPECT_EQ(occmat.data()[iat][l][1].nc, 5);

    EXPECT_EQ(occmat.data_save()[iat][l].size(), 2u);
    EXPECT_EQ(occmat.data_save()[iat][l][0].nr, 5);
    EXPECT_EQ(occmat.data_save()[iat][l][0].nc, 5);
    EXPECT_EQ(occmat.data_save()[iat][l][1].nr, 5);
    EXPECT_EQ(occmat.data_save()[iat][l][1].nc, 5);
}

// ----------------------------------------------------------------------
// init: nspin=2 layout
// ----------------------------------------------------------------------
TEST_F(OccMatrixTest, InitCollinearNspin2)
{
    make_cell({2}, {2}, {{0, 0, 1}}, 1);
    OccupationMatrix occmat;
    occmat.init(cell, {2}, 2, 1);

    EXPECT_EQ(occmat.nspin(), 2);
    EXPECT_EQ(occmat.npol(), 1);

    const int iat = 0;
    const int l = 2;
    EXPECT_EQ(occmat.data()[iat][l].size(), 2u);
    EXPECT_EQ(occmat.data()[iat][l][0].nr, 5);
    EXPECT_EQ(occmat.data()[iat][l][0].nc, 5);
    EXPECT_EQ(occmat.data()[iat][l][1].nr, 5);
    EXPECT_EQ(occmat.data()[iat][l][1].nc, 5);

    EXPECT_EQ(occmat.data_save()[iat][l].size(), 2u);
    EXPECT_EQ(occmat.data_save()[iat][l][0].nr, 5);
    EXPECT_EQ(occmat.data_save()[iat][l][0].nc, 5);
    EXPECT_EQ(occmat.data_save()[iat][l][1].nr, 5);
    EXPECT_EQ(occmat.data_save()[iat][l][1].nc, 5);
}

// ----------------------------------------------------------------------
// init: nspin=4 layout
// ----------------------------------------------------------------------
TEST_F(OccMatrixTest, InitNoncollinearLayout)
{
    make_cell({2}, {2}, {{0, 0, 1}}, 2);
    OccupationMatrix occmat;
    occmat.init(cell, {2}, 4, 2);

    EXPECT_EQ(occmat.nspin(), 4);
    EXPECT_EQ(occmat.npol(), 2);

    const int iat = 0;
    const int l = 2;
    EXPECT_EQ(occmat.data()[iat][l].size(), 1u);
    EXPECT_EQ(occmat.data()[iat][l][0].nr, 10);
    EXPECT_EQ(occmat.data()[iat][l][0].nc, 10);

    EXPECT_EQ(occmat.data_save()[iat][l].size(), 1u);
    EXPECT_EQ(occmat.data_save()[iat][l][0].nr, 10);
    EXPECT_EQ(occmat.data_save()[iat][l][0].nc, 10);
}

// ----------------------------------------------------------------------
// init: uncorrelated type is skipped and does not pollute the rest
// ----------------------------------------------------------------------
TEST_F(OccMatrixTest, InitSkipsUncorrelatedType)
{
    make_cell({2, -1}, {2, 2}, {{0, 0, 1}, {0, 0, 1}}, 1);
    OccupationMatrix occmat;
    occmat.init(cell, {2, -1}, 1, 1);

    // correlated type 0: allocated
    EXPECT_EQ(occmat.data()[0][2].size(), 2u);
    EXPECT_EQ(occmat.data()[0][2][0].nr, 5);
    EXPECT_EQ(occmat.data_save()[0][2].size(), 2u);

    // uncorrelated type 1: no valid occupation blocks
    EXPECT_EQ(occmat.data()[1][2].size(), 0u);
    EXPECT_EQ(occmat.data_save()[1][2].size(), 0u);
}

// ----------------------------------------------------------------------
// init: corr_iwt points to the first radial channel (n=0)
// ----------------------------------------------------------------------
TEST_F(OccMatrixTest, InitBuildsCorrIwtFromFirstRadialChannel)
{
    // one type, one atom, l=2 with two radial channels (n=0 and n=1)
    make_cell({2}, {2}, {{0, 0, 2}}, 1);
    OccupationMatrix occmat;
    occmat.init(cell, {2}, 1, 1);

    const int iat = 0;
    const int l = 2;
    const int m_size = 5;

    // global indices: atom 0 starts at 0, nw = 2*5 = 10
    // n=0: iw = 0..4, n=1: iw = 5..9
    for (int m = 0; m < m_size; m++)
    {
        EXPECT_EQ(occmat.corr_iwt(iat, l, m, 0), m);
    }
}

TEST_F(OccMatrixTest, InitBuildsCorrIwtForLaterNoncollinearAtom)
{
    // Atom 0 contributes one s orbital with two polarizations, so atom 1
    // starts at global index 2. Its d orbitals are interleaved by ipol.
    make_cell({-1, 2}, {0, 2}, {{1}, {0, 0, 1}}, 2);
    OccupationMatrix occmat;
    occmat.init(cell, {-1, 2}, 4, 2);

    const int iat = 1;
    const int l = 2;
    for (int m = 0; m < 5; m++)
    {
        EXPECT_EQ(occmat.corr_iwt(iat, l, m, 0), 2 + 2 * m);
        EXPECT_EQ(occmat.corr_iwt(iat, l, m, 1), 2 + 2 * m + 1);
    }
}

// ----------------------------------------------------------------------
// write_save_to_flat: nspin=1
// ----------------------------------------------------------------------
TEST_F(OccMatrixTest, WriteSaveToFlatNspin1)
{
    make_cell({2}, {2}, {{0, 0, 1}}, 1);
    OccupationMatrix occmat;
    occmat.init(cell, {2}, 1, 1);

    const int block = 25;
    std::vector<int> index = {0};
    std::vector<double> uom_save(block, -1.0);

    for (int m = 0; m < block; m++)
    {
        occmat.data_save()[0][2][0].c[m] = 10.0 + m;
    }

    occmat.write_save_to_flat(cell, {2}, index, uom_save);

    for (int m = 0; m < block; m++)
    {
        EXPECT_DOUBLE_EQ(uom_save[index[0] + m], 10.0 + m);
    }
}

// ----------------------------------------------------------------------
// write_save_to_flat: nspin=2 split layout
// ----------------------------------------------------------------------
TEST_F(OccMatrixTest, WriteSaveToFlatNspin2SplitLayout)
{
    // two atoms, both correlated d
    make_cell({2, 2}, {2, 2}, {{0, 0, 1}, {0, 0, 1}}, 1);
    OccupationMatrix occmat;
    occmat.init(cell, {2, 2}, 2, 1);

    const int block = 25;
    const int half = cell.nat * block;
    std::vector<int> index = {0, block};
    std::vector<double> uom_save(2 * half, -1.0);

    // four distinct data sets
    for (int m = 0; m < block; m++)
    {
        occmat.data_save()[0][2][0].c[m] = 1000.0 + m; // atom0 up
        occmat.data_save()[1][2][0].c[m] = 2000.0 + m; // atom1 up
        occmat.data_save()[0][2][1].c[m] = 3000.0 + m; // atom0 down
        occmat.data_save()[1][2][1].c[m] = 4000.0 + m; // atom1 down
    }

    occmat.write_save_to_flat(cell, {2, 2}, index, uom_save);

    for (int m = 0; m < block; m++)
    {
        EXPECT_DOUBLE_EQ(uom_save[0 * block + m], 1000.0 + m);
        EXPECT_DOUBLE_EQ(uom_save[1 * block + m], 2000.0 + m);
        EXPECT_DOUBLE_EQ(uom_save[half + 0 * block + m], 3000.0 + m);
        EXPECT_DOUBLE_EQ(uom_save[half + 1 * block + m], 4000.0 + m);
    }
}

// ----------------------------------------------------------------------
// write_save_to_flat: nspin=4 full Pauli block
// ----------------------------------------------------------------------
TEST_F(OccMatrixTest, WriteSaveToFlatNspin4)
{
    make_cell({2}, {2}, {{0, 0, 1}}, 2);
    OccupationMatrix occmat;
    occmat.init(cell, {2}, 4, 2);

    const int block = 100; // 10 x 10
    std::vector<int> index = {0};
    std::vector<double> uom_save(block, -1.0);

    for (int m = 0; m < block; m++)
    {
        occmat.data_save()[0][2][0].c[m] = 500.0 + m;
    }

    occmat.write_save_to_flat(cell, {2}, index, uom_save);

    for (int m = 0; m < block; m++)
    {
        EXPECT_DOUBLE_EQ(uom_save[m], 500.0 + m);
    }
}

// ----------------------------------------------------------------------
// write_save_to_flat: empty buffer is a no-op
// ----------------------------------------------------------------------
TEST_F(OccMatrixTest, WriteSaveToFlatEmptyBufferIsNoOp)
{
    make_cell({2}, {2}, {{0, 0, 1}}, 1);
    OccupationMatrix occmat;
    occmat.init(cell, {2}, 1, 1);

    std::vector<int> index = {0};
    std::vector<double> uom_save;

    occmat.write_save_to_flat(cell, {2}, index, uom_save);
    EXPECT_TRUE(uom_save.empty());
}
