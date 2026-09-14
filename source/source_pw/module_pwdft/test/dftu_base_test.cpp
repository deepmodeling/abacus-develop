/**********************************************
 *  Unit tests for Plus_U_Base::init_base.
 *
 *  Focus: the Yukawa state must follow the
 *  yukawa_potential argument on every call,
 *  including a true -> false re-initialization
 *  (before_scf() -> setup_pot() may call
 *  init_base() repeatedly on the same object).
 ***********************************************/

#include "source_pw/module_pwdft/dftu_base.h"

#include "source_cell/atom_spec.h"
#include "source_cell/unitcell.h"

#include "gtest/gtest.h"

#include <vector>
#include <numeric>

class DFTUBaseTest : public testing::Test
{
  protected:
    UnitCell ucell;

    void SetUp() override
    {
        // Minimal one-atom cell: d channel available (nwl = 2),
        // one chi for each of s / p / d, so nw = 1 + 3 + 5 = 9.
        const int nw = 9;

        ucell.ntype = 1;
        ucell.nat = 1;
        ucell.atoms = new Atom[ucell.ntype];
        ucell.iat2it = new int[ucell.nat];
        ucell.iat2ia = new int[ucell.nat];
        ucell.atoms[0].tau.resize(ucell.nat);
        ucell.atoms[0].taud.resize(ucell.nat);
        ucell.itia2iat.create(ucell.ntype, ucell.nat);
        for (int iat = 0; iat < ucell.nat; iat++)
        {
            ucell.iat2it[iat] = 0;
            ucell.iat2ia[iat] = iat;
            ucell.itia2iat(0, iat) = iat;
            ucell.atoms[0].tau[iat] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
            ucell.atoms[0].taud[iat] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
        }
        ucell.atoms[0].na = 1;
        ucell.atoms[0].label = "Fe";
        ucell.atoms[0].nwl = 2;
        ucell.atoms[0].l_nchi = {1, 1, 1};
        ucell.atoms[0].nw = nw;
        ucell.atoms[0].iw2l.resize(nw);
        ucell.atoms[0].iw2n.resize(nw);
        ucell.atoms[0].iw2m.resize(nw);
        int iw = 0;
        for (int l = 0; l <= ucell.atoms[0].nwl; l++)
        {
            for (int m = 0; m < 2 * l + 1; m++)
            {
                ucell.atoms[0].iw2l[iw] = l;
                ucell.atoms[0].iw2n[iw] = 0;
                ucell.atoms[0].iw2m[iw] = m;
                iw++;
            }
        }
        ucell.set_iat2iwt(1);
    }

    void TearDown() override
    {
        // set_atom_flag is false, so ~UnitCell() skips atoms but frees
        // iat2it / iat2ia itself; only atoms must be deleted here.
        delete[] ucell.atoms;
    }

    /// Call init_base with the given Yukawa switch on a fresh d orbital
    void init_dftu(Plus_U_Base& dftu, const bool yukawa_potential)
    {
        const std::vector<int> l_channel = {2};
        const std::vector<double> hubbard_u = {0.0};
        dftu.init_base(ucell,
                       1,                      // npol
                       2,                      // nspin
                       l_channel,
                       yukawa_potential,
                       0.5,                    // yukawa_lambda
                       "",                     // global_readin_dir
                       "",                     // global_out_dir
                       "none",                 // init_chg
                       "cpu",                  // device
                       hubbard_u,
                       0.0,                    // uramping
                       0,                      // occ_mat_ctrl
                       0);                     // mixing_dftu
    }
};

/// After a true -> false re-initialization the Yukawa object must be
/// released so that use_yukawa() reflects the latest argument.
TEST_F(DFTUBaseTest, InitBaseYukawaTrueThenFalseClearsState)
{
    Plus_U_Base dftu;

    init_dftu(dftu, true);
    EXPECT_TRUE(dftu.use_yukawa());

    init_dftu(dftu, false);
    EXPECT_FALSE(dftu.use_yukawa());
}

/// A false -> true re-initialization must create the Yukawa object.
TEST_F(DFTUBaseTest, InitBaseYukawaFalseThenTrueCreatesObject)
{
    Plus_U_Base dftu;

    init_dftu(dftu, false);
    EXPECT_FALSE(dftu.use_yukawa());

    init_dftu(dftu, true);
    EXPECT_TRUE(dftu.use_yukawa());
}

// =====================================================================
// uterm_mat_index calculation
//
// nspin=1: offset = sum(tlp1^2), total = sum(all tlp1^2)
// nspin=2: same per-spin-channel, then pot_index *= 2 (split layout)
// nspin=4: offset = sum((tlp1*npol)^2), each atom = 4*tlp1^2
// =====================================================================

class EffPotIndexTest : public ::testing::Test
{
  protected:
    struct AtomSpec { int l; int na; }; // correlated orbital l, number of atoms
    std::vector<int> uterm_mat_index;
    int pot_index;

    void compute_indices(const std::vector<AtomSpec>& atoms, int nspin)
    {
        pot_index = 0;
        uterm_mat_index.resize(atoms.size());

        for (size_t i = 0; i < atoms.size(); i++)
        {
            int tlp1 = 2 * atoms[i].l + 1;
            int tlp1_npol = tlp1 * (nspin == 4 ? 2 : 1);

            if (nspin == 4)
            {
                uterm_mat_index[i] = pot_index;
                pot_index += tlp1_npol * tlp1_npol;
            }
            else
            {
                uterm_mat_index[i] = pot_index;
                pot_index += tlp1 * tlp1;
            }
        }

        if (nspin == 2)
            pot_index *= 2;
    }
};

TEST_F(EffPotIndexTest, Nspin1_MixedOrbitals)
{
    // 3 atoms: p(l=1), d(l=2), p(l=1)
    std::vector<AtomSpec> atoms = {{1, 1}, {2, 1}, {1, 1}};
    compute_indices(atoms, 1);

    // p: 9, d: 25, p: 9
    EXPECT_EQ(uterm_mat_index[0], 0);
    EXPECT_EQ(uterm_mat_index[1], 9);
    EXPECT_EQ(uterm_mat_index[2], 34);
    EXPECT_EQ(pot_index, 43); // 9 + 25 + 9
}

TEST_F(EffPotIndexTest, Nspin2and4_SplitAndPauli)
{
    // nspin=2: 2 d-atoms, split layout [up | dn]
    std::vector<AtomSpec> atoms2 = {{2, 1}, {2, 1}};
    compute_indices(atoms2, 2);
    EXPECT_EQ(uterm_mat_index[0], 0);
    EXPECT_EQ(uterm_mat_index[1], 25);
    EXPECT_EQ(pot_index, 100); // (25 + 25) * 2

    // nspin=4: d + p atoms, Pauli blocks
    std::vector<AtomSpec> atoms4 = {{2, 1}, {1, 1}};
    compute_indices(atoms4, 4);
    EXPECT_EQ(uterm_mat_index[0], 0);    // d: (5*2)^2 = 100
    EXPECT_EQ(uterm_mat_index[1], 100);  // p: (3*2)^2 = 36
    EXPECT_EQ(pot_index, 136);
}

// =====================================================================
// copy_occ_mat <-> set_occ_mat roundtrip
//
// Tests the bidirectional conversion between nested occ_mat matrix
// and flat uom_array/uom_save arrays for all 3 nspin modes.
// =====================================================================

struct Matrix2D {
    int nr, nc;
    std::vector<double> data;
    Matrix2D() : nr(0), nc(0), data() {}
    Matrix2D(int r, int c) : nr(r), nc(c), data(r * c, 0.0) {}
    double& operator()(int i, int j) { return data[i * nc + j]; }
    const double& operator()(int i, int j) const { return data[i * nc + j]; }
};

static void copy_occ_mat_to_flat(
    const std::vector<Matrix2D>& occ_mat_up,
    const std::vector<Matrix2D>& occ_mat_dn,
    std::vector<double>& uom_save,
    const std::vector<int>& uterm_mat_index,
    int nspin)
{
    if (nspin == 4)
    {
        for (size_t iat = 0; iat < occ_mat_up.size(); iat++)
        {
            int size = occ_mat_up[iat].nr * occ_mat_up[iat].nc;
            for (int mm = 0; mm < size; mm++)
                uom_save[uterm_mat_index[iat] + mm] = occ_mat_up[iat].data[mm];
        }
    }
    else if (nspin == 2) // split layout: [up | dn]
    {
        int half_size = uom_save.size() / 2;
        for (size_t iat = 0; iat < occ_mat_up.size(); iat++)
        {
            int size = occ_mat_up[iat].nr * occ_mat_up[iat].nc;
            for (int mm = 0; mm < size; mm++)
            {
                uom_save[uterm_mat_index[iat] + mm] = occ_mat_up[iat].data[mm];
                uom_save[half_size + uterm_mat_index[iat] + mm] = occ_mat_dn[iat].data[mm];
            }
        }
    }
    else // nspin=1: single spin channel
    {
        for (size_t iat = 0; iat < occ_mat_up.size(); iat++)
        {
            int size = occ_mat_up[iat].nr * occ_mat_up[iat].nc;
            for (int mm = 0; mm < size; mm++)
                uom_save[uterm_mat_index[iat] + mm] = occ_mat_up[iat].data[mm];
        }
    }
}

static void set_occ_mat_from_flat(
    const std::vector<double>& uom_array,
    std::vector<Matrix2D>& occ_mat_up,
    std::vector<Matrix2D>& occ_mat_dn,
    const std::vector<int>& uterm_mat_index,
    int nspin)
{
    if (nspin == 4)
    {
        for (size_t iat = 0; iat < occ_mat_up.size(); iat++)
        {
            int size = occ_mat_up[iat].nr * occ_mat_up[iat].nc;
            for (int mm = 0; mm < size; mm++)
                occ_mat_up[iat].data[mm] = uom_array[uterm_mat_index[iat] + mm];
        }
    }
    else if (nspin == 2)
    {
        int half_size = uom_array.size() / 2;
        for (size_t iat = 0; iat < occ_mat_up.size(); iat++)
        {
            int size = occ_mat_up[iat].nr * occ_mat_up[iat].nc;
            for (int mm = 0; mm < size; mm++)
            {
                occ_mat_up[iat].data[mm] = uom_array[uterm_mat_index[iat] + mm];
                occ_mat_dn[iat].data[mm] = uom_array[half_size + uterm_mat_index[iat] + mm];
            }
        }
    }
    else // nspin=1
    {
        for (size_t iat = 0; iat < occ_mat_up.size(); iat++)
        {
            int size = occ_mat_up[iat].nr * occ_mat_up[iat].nc;
            for (int mm = 0; mm < size; mm++)
                occ_mat_up[iat].data[mm] = uom_array[uterm_mat_index[iat] + mm];
        }
    }
}

class OccMatRoundtripTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
};

TEST_F(OccMatRoundtripTest, Nspin1and2_SingleAndSplitLayout)
{
    // nspin=1: single atom d-orbital roundtrip
    const int l = 2;
    const int size = (2 * l + 1) * (2 * l + 1); // 25

    std::vector<Matrix2D> occ_mat_up(1, Matrix2D(2 * l + 1, 2 * l + 1));
    std::vector<Matrix2D> occ_mat_dn(1, Matrix2D(2 * l + 1, 2 * l + 1));
    for (int i = 0; i < size; i++)
        occ_mat_up[0].data[i] = static_cast<double>(i + 1);

    std::vector<int> uterm_mat_index = {0};
    std::vector<double> uom_save(size, 0.0);
    copy_occ_mat_to_flat(occ_mat_up, occ_mat_dn, uom_save, uterm_mat_index, 1);
    set_occ_mat_from_flat(uom_save, occ_mat_up, occ_mat_dn, uterm_mat_index, 1);
    for (int i = 0; i < size; i++)
        EXPECT_DOUBLE_EQ(occ_mat_up[0].data[i], static_cast<double>(i + 1));

    // nspin=2: split layout [up | dn] with distinct values
    const int total = size * 2;
    for (int i = 0; i < size; i++)
    {
        occ_mat_up[0].data[i] = static_cast<double>(i + 1);
        occ_mat_dn[0].data[i] = static_cast<double>(i + 100);
    }
    uom_save.assign(total, 0.0);
    copy_occ_mat_to_flat(occ_mat_up, occ_mat_dn, uom_save, uterm_mat_index, 2);
    // Verify split layout
    for (int i = 0; i < size; i++)
    {
        EXPECT_DOUBLE_EQ(uom_save[i], static_cast<double>(i + 1));
        EXPECT_DOUBLE_EQ(uom_save[size + i], static_cast<double>(i + 100));
    }
    set_occ_mat_from_flat(uom_save, occ_mat_up, occ_mat_dn, uterm_mat_index, 2);
    for (int i = 0; i < size; i++)
    {
        EXPECT_DOUBLE_EQ(occ_mat_up[0].data[i], static_cast<double>(i + 1));
        EXPECT_DOUBLE_EQ(occ_mat_dn[0].data[i], static_cast<double>(i + 100));
    }
}

TEST_F(OccMatRoundtripTest, Nspin4_PauliBlocks)
{
    // 2 atoms: d(l=2), p(l=1)
    struct AtomSpec { int l; };
    std::vector<AtomSpec> specs = {{2}, {1}};
    int npol = 2;

    std::vector<int> sizes;
    for (auto& s : specs)
    {
        int tlp1 = 2 * s.l + 1;
        sizes.push_back((tlp1 * npol) * (tlp1 * npol));
    }
    int total = std::accumulate(sizes.begin(), sizes.end(), 0);

    std::vector<int> uterm_mat_index(specs.size());
    int offset = 0;
    for (size_t i = 0; i < specs.size(); i++)
    {
        uterm_mat_index[i] = offset;
        offset += sizes[i];
    }

    std::vector<Matrix2D> occ_mat(specs.size());
    for (size_t i = 0; i < specs.size(); i++)
    {
        int dim = (2 * specs[i].l + 1) * npol;
        occ_mat[i] = Matrix2D(dim, dim);
        for (int j = 0; j < sizes[i]; j++)
            occ_mat[i].data[j] = static_cast<double>(i * 1000 + j + 1);
    }

    std::vector<double> uom_array(total, 0.0);
    std::vector<Matrix2D> occ_mat_dn(specs.size()); // unused for nspin=4

    copy_occ_mat_to_flat(occ_mat, occ_mat_dn, uom_array, uterm_mat_index, 4);
    set_occ_mat_from_flat(uom_array, occ_mat, occ_mat_dn, uterm_mat_index, 4);

    for (size_t i = 0; i < specs.size(); i++)
        for (int j = 0; j < sizes[i]; j++)
            EXPECT_DOUBLE_EQ(occ_mat[i].data[j], static_cast<double>(i * 1000 + j + 1));
}
