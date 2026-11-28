#include "source_estate/occ_mixer.h"

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
 * Unit tests for OccMatMixer.
 *
 * Covered:
 * 1. mix_plain: occ = beta*occ + (1-beta)*occ_save on nested blocks
 *    (nspin=1 and nspin=2, including the uncorrelated-atom skip).
 * 2. flat buffer roundtrip: collect -> write_back reproduces the nested
 *    occupation matrix for the split (nspin=2) and spinor (nspin=4) layouts.
 *
 * A minimal UnitCell with one type, two atoms, a single correlated
 * d-orbital (l=2, one radial channel) is built in SetUp.
 ***********************************************************************/

class OccMatMixerTest : public ::testing::Test
{
  protected:
    // one correlated d-orbital per atom; l=2 -> 5x5 block per spin channel
    static const int l_corr = 2;
    static const int m_size = 2 * l_corr + 1; // 5
    static const int block = m_size * m_size; // 25

    void SetUp() override
    {
        // minimal unit cell: 1 type, 2 atoms, both correlated d-atoms
        cell.ntype = 1;
        cell.nat = 2;
        atoms_storage.resize(cell.ntype);
        cell.atoms = atoms_storage.data();
        cell.iat2it.resize(cell.nat);
        cell.iat2ia.resize(cell.nat);
        cell.itia2iat.create(cell.ntype, cell.nat);
        cell.atoms[0].na = 2;
        cell.atoms[0].nwl = 2;      // max angular momentum present
        cell.atoms[0].l_nchi.resize(cell.atoms[0].nwl + 1);
        cell.atoms[0].l_nchi[0] = 0;
        cell.atoms[0].l_nchi[1] = 0;
        cell.atoms[0].l_nchi[2] = 1; // one d radial channel
        for (int iat = 0; iat < cell.nat; iat++)
        {
            cell.iat2it[iat] = 0;
            cell.iat2ia[iat] = iat;
            cell.itia2iat(0, iat) = iat;
        }

        l_channel = {l_corr};

        // per-atom offset table matching Plus_U_Base::init_base layout
        flat_index.resize(cell.nat);
        flat_index[0] = 0;
        flat_index[1] = block;
    }

    // total flat size: nspin=2 doubles (split [up | dn]); nspin=1/4 single
    int total_size(const int nspin) const
    {
        const int per_spin = cell.nat * block;
        return (nspin == 2) ? 2 * per_spin : per_spin;
    }

    UnitCell cell;
    std::vector<Atom> atoms_storage;  ///< backing store for cell.atoms
    std::vector<int> l_channel;
    std::vector<int> flat_index;
};

const int OccMatMixerTest::l_corr;
const int OccMatMixerTest::m_size;
const int OccMatMixerTest::block;

// ----------------------------------------------------------------------
// mix_plain: nested-matrix linear mixing
// ----------------------------------------------------------------------
TEST_F(OccMatMixerTest, MixPlainNspin1)
{
    const int nspin = 1;
    const int npol = 1;
    OccupationMatrix occmat;
    occmat.init(cell, l_channel, nspin, npol);

    OccMatMixer mixer;
    mixer.init(&cell, &l_channel, &flat_index, nspin, total_size(nspin));

    // fill occ and occ_save with known distinct values
    for (int iat = 0; iat < cell.nat; iat++)
    {
        for (int m = 0; m < block; m++)
        {
            occmat.data()[iat][l_corr][0].c[m] = 1.0 + m;
            occmat.data_save()[iat][l_corr][0].c[m] = 100.0 + m;
        }
    }

    const double beta = 0.25;
    mixer.mix_plain(occmat, beta);

    for (int iat = 0; iat < cell.nat; iat++)
    {
        for (int m = 0; m < block; m++)
        {
            const double expect = (1.0 + m) * beta + (100.0 + m) * (1.0 - beta);
            EXPECT_DOUBLE_EQ(occmat.data()[iat][l_corr][0].c[m], expect);
        }
    }
}

TEST_F(OccMatMixerTest, MixPlainNspin2BothChannels)
{
    const int nspin = 2;
    const int npol = 1;
    OccupationMatrix occmat;
    occmat.init(cell, l_channel, nspin, npol);

    OccMatMixer mixer;
    mixer.init(&cell, &l_channel, &flat_index, nspin, total_size(nspin));

    for (int iat = 0; iat < cell.nat; iat++)
    {
        for (int is = 0; is < 2; is++)
        {
            for (int m = 0; m < block; m++)
            {
                occmat.data()[iat][l_corr][is].c[m] = 2.0 + is + m;
                occmat.data_save()[iat][l_corr][is].c[m] = 50.0 + is + m;
            }
        }
    }

    const double beta = 0.5;
    mixer.mix_plain(occmat, beta);

    for (int iat = 0; iat < cell.nat; iat++)
    {
        for (int is = 0; is < 2; is++)
        {
            for (int m = 0; m < block; m++)
            {
                const double expect = (2.0 + is + m) * beta + (50.0 + is + m) * (1.0 - beta);
                EXPECT_DOUBLE_EQ(occmat.data()[iat][l_corr][is].c[m], expect);
            }
        }
    }
}

// ----------------------------------------------------------------------
// flat buffer roundtrip: collect then write_back must reproduce occ
// ----------------------------------------------------------------------
TEST_F(OccMatMixerTest, FlatRoundtripNspin2)
{
    const int nspin = 2;
    const int npol = 1;
    OccupationMatrix occmat;
    occmat.init(cell, l_channel, nspin, npol);

    OccMatMixer mixer;
    mixer.init(&cell, &l_channel, &flat_index, nspin, total_size(nspin));
    EXPECT_EQ(mixer.flat_size(), total_size(nspin));

    // distinct value per (iat, spin, m) to detect layout mistakes
    for (int iat = 0; iat < cell.nat; iat++)
    {
        for (int is = 0; is < 2; is++)
        {
            for (int m = 0; m < block; m++)
            {
                occmat.data()[iat][l_corr][is].c[m] =
                    1000.0 * iat + 100.0 * is + m;
            }
        }
    }

    mixer.collect(occmat);   // occ -> uom_
    const int half_size = cell.nat * block;
    for (int iat = 0; iat < cell.nat; ++iat)
    {
        for (int m = 0; m < block; ++m)
        {
            const int up_index = flat_index[iat] + m;
            const int down_index = half_size + up_index;
            EXPECT_DOUBLE_EQ(mixer.uom()[up_index], 1000.0 * iat + m);
            EXPECT_DOUBLE_EQ(mixer.uom()[down_index], 1000.0 * iat + 100.0 + m);
        }
    }
    // scramble the nested matrix, then restore it from the flat buffer
    occmat.zero(cell, l_channel);
    mixer.write_back(occmat); // uom_ -> occ

    for (int iat = 0; iat < cell.nat; iat++)
    {
        for (int is = 0; is < 2; is++)
        {
            for (int m = 0; m < block; m++)
            {
                EXPECT_DOUBLE_EQ(occmat.data()[iat][l_corr][is].c[m],
                                 1000.0 * iat + 100.0 * is + m);
            }
        }
    }
}

TEST_F(OccMatMixerTest, FlatRoundtripNspin4AllPauliBlocks)
{
    const int nspin = 4;
    const int npol = 2;
    const int spinor_block = block * npol * npol;
    const int flat_size = cell.nat * spinor_block;
    flat_index[1] = spinor_block;

    OccupationMatrix occmat;
    occmat.init(cell, l_channel, nspin, npol);

    OccMatMixer mixer;
    mixer.init(&cell, &l_channel, &flat_index, nspin, flat_size);

    for (int iat = 0; iat < cell.nat; iat++)
    {
        for (int m = 0; m < spinor_block; m++)
        {
            occmat.mat(iat, l_corr, 0).c[m] = 1000.0 * iat + m + 1.0;
        }
    }

    mixer.collect(occmat);
    for (int iat = 0; iat < cell.nat; iat++)
    {
        for (int m = 0; m < spinor_block; m++)
        {
            EXPECT_DOUBLE_EQ(mixer.uom()[flat_index[iat] + m],
                             1000.0 * iat + m + 1.0);
        }
    }

    occmat.zero(cell, l_channel);
    mixer.write_back(occmat);
    for (int iat = 0; iat < cell.nat; iat++)
    {
        for (int m = 0; m < spinor_block; m++)
        {
            EXPECT_DOUBLE_EQ(occmat.mat(iat, l_corr, 0).c[m],
                             1000.0 * iat + m + 1.0);
        }
    }
}

// ----------------------------------------------------------------------
// begin_iter flattens the saved (not the live) occupation matrix
// ----------------------------------------------------------------------
TEST_F(OccMatMixerTest, BeginIterFlattensSave)
{
    const int nspin = 1;
    const int npol = 1;
    OccupationMatrix occmat;
    occmat.init(cell, l_channel, nspin, npol);

    OccMatMixer mixer;
    mixer.init(&cell, &l_channel, &flat_index, nspin, total_size(nspin));

    for (int m = 0; m < block; m++)
    {
        occmat.data_save()[0][l_corr][0].c[m] = 7.0 + m;
    }

    mixer.begin_iter(occmat);

    for (int m = 0; m < block; m++)
    {
        EXPECT_DOUBLE_EQ(mixer.uom_save()[flat_index[0] + m], 7.0 + m);
    }
}

TEST_F(OccMatMixerTest, OnsiteFirstProjectionAndIndependentSpinHistory)
{
    OccupationMatrix occupation;
    occupation.init(cell, l_channel, 2, 1);
    OccMatMixer mixer;
    const int size = total_size(2);
    mixer.init(&cell, &l_channel, &flat_index, 2, size);
    std::vector<double> up(block, 0.2);
    std::vector<double> down(block, 0.6);
    mixer.mix_onsite_block(0, 0, up, 0.25);
    mixer.mix_onsite_block(0, 1, down, 0.25);
    EXPECT_DOUBLE_EQ(up[0], 0.2);
    EXPECT_DOUBLE_EQ(down[0], 0.6);
    EXPECT_TRUE(mixer.mixed_block(0, 0, block).empty());
    up.assign(block, 0.8);
    down.assign(block, 0.2);
    mixer.mix_onsite_block(0, 0, up, 0.25);
    mixer.mix_onsite_block(0, 1, down, 0.25);
    EXPECT_NEAR(up[0], 0.35, 1e-15);
    EXPECT_NEAR(down[0], 0.5, 1e-15);
    EXPECT_DOUBLE_EQ(mixer.raw_block(0, 0, block)[0], 0.8);
    EXPECT_EQ(mixer.mixed_block(0, 0, block), up);
    // Another atom has no previous matrix, and must not inherit atom 0's history.
    std::vector<double> other(block, 0.9);
    mixer.mix_onsite_block(1, 0, other, 0.25);
    EXPECT_DOUBLE_EQ(other[0], 0.9);
}

TEST_F(OccMatMixerTest, SeededOnsiteMatrixAndPwRawSnapshot)
{
    OccupationMatrix occupation;
    occupation.init(cell, l_channel, 2, 1);
    occupation.zero(cell, l_channel);
    std::vector<double> input(block, 0.4);
    occupation.set_flat(0, l_corr, 0, input);
    OccMatMixer mixer;
    const int size = total_size(2);
    mixer.init(&cell, &l_channel, &flat_index, 2, size);
    mixer.seed_save(occupation);
    std::vector<double> raw(block, 0.8);
    mixer.mix_onsite_block(0, 0, raw, 0.5);
    EXPECT_NEAR(raw[0], 0.6, 1e-15);
    occupation.copy_to_save(cell, l_channel);
    mixer.begin_iter(occupation);
    occupation.set_flat(0, l_corr, 0, raw);
    mixer.collect(occupation);
    mixer.uom()[0] = 0.1;
    mixer.record_mixing(true);
    mixer.write_back(occupation);
    EXPECT_NEAR(mixer.raw_block(0, 0, block)[0], 0.6, 1e-15);
    EXPECT_DOUBLE_EQ(mixer.mixed_block(0, 0, block)[0], 0.1);
}
