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
 *    occupation matrix for the split (nspin=2) layout.
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
        // UnitCell exposes raw pointers; Statistics::~Statistics() owns and
        // delete[]s iat2it/iat2ia, so they must come from new[] here.
        cell.atoms = atoms_storage.data();
        cell.iat2it = new int[cell.nat];
        cell.iat2ia = new int[cell.nat];
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

        orbital_corr = {l_corr};

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
    std::vector<int> orbital_corr;
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
    occmat.init(cell, orbital_corr, nspin, npol);

    OccMatMixer mixer;
    mixer.init(&cell, &orbital_corr, &flat_index, nspin, total_size(nspin));

    // fill occ and occ_save with known distinct values
    for (int iat = 0; iat < cell.nat; iat++)
    {
        for (int m = 0; m < block; m++)
        {
            occmat.data()[iat][l_corr][0][0].c[m] = 1.0 + m;
            occmat.data_save()[iat][l_corr][0][0].c[m] = 100.0 + m;
        }
    }

    const double beta = 0.25;
    mixer.mix_plain(occmat, beta);

    for (int iat = 0; iat < cell.nat; iat++)
    {
        for (int m = 0; m < block; m++)
        {
            const double expect = (1.0 + m) * beta + (100.0 + m) * (1.0 - beta);
            EXPECT_DOUBLE_EQ(occmat.data()[iat][l_corr][0][0].c[m], expect);
        }
    }
}

TEST_F(OccMatMixerTest, MixPlainNspin2BothChannels)
{
    const int nspin = 2;
    const int npol = 1;
    OccupationMatrix occmat;
    occmat.init(cell, orbital_corr, nspin, npol);

    OccMatMixer mixer;
    mixer.init(&cell, &orbital_corr, &flat_index, nspin, total_size(nspin));

    for (int iat = 0; iat < cell.nat; iat++)
    {
        for (int is = 0; is < 2; is++)
        {
            for (int m = 0; m < block; m++)
            {
                occmat.data()[iat][l_corr][0][is].c[m] = 2.0 + is + m;
                occmat.data_save()[iat][l_corr][0][is].c[m] = 50.0 + is + m;
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
                EXPECT_DOUBLE_EQ(occmat.data()[iat][l_corr][0][is].c[m], expect);
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
    occmat.init(cell, orbital_corr, nspin, npol);

    OccMatMixer mixer;
    mixer.init(&cell, &orbital_corr, &flat_index, nspin, total_size(nspin));
    EXPECT_EQ(mixer.flat_size(), total_size(nspin));

    // distinct value per (iat, spin, m) to detect layout mistakes
    for (int iat = 0; iat < cell.nat; iat++)
    {
        for (int is = 0; is < 2; is++)
        {
            for (int m = 0; m < block; m++)
            {
                occmat.data()[iat][l_corr][0][is].c[m] =
                    1000.0 * iat + 100.0 * is + m;
            }
        }
    }

    mixer.collect(occmat);   // occ -> uom_
    // scramble the nested matrix, then restore it from the flat buffer
    occmat.zero(cell, orbital_corr);
    mixer.write_back(occmat); // uom_ -> occ

    for (int iat = 0; iat < cell.nat; iat++)
    {
        for (int is = 0; is < 2; is++)
        {
            for (int m = 0; m < block; m++)
            {
                EXPECT_DOUBLE_EQ(occmat.data()[iat][l_corr][0][is].c[m],
                                 1000.0 * iat + 100.0 * is + m);
            }
        }
    }
}

// ----------------------------------------------------------------------
// begin_iter/seed_save flatten the saved (not the live) occupation matrix
// ----------------------------------------------------------------------
TEST_F(OccMatMixerTest, BeginIterFlattensSave)
{
    const int nspin = 1;
    const int npol = 1;
    OccupationMatrix occmat;
    occmat.init(cell, orbital_corr, nspin, npol);

    OccMatMixer mixer;
    mixer.init(&cell, &orbital_corr, &flat_index, nspin, total_size(nspin));

    for (int m = 0; m < block; m++)
    {
        occmat.data_save()[0][l_corr][0][0].c[m] = 7.0 + m;
    }

    mixer.begin_iter(occmat);

    for (int m = 0; m < block; m++)
    {
        EXPECT_DOUBLE_EQ(mixer.uom_save()[flat_index[0] + m], 7.0 + m);
    }
}
