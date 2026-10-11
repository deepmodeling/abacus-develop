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
#include "source_pw/module_pwdft/dftu_base_io.h"

#include "source_cell/atom_spec.h"
#include "source_cell/unitcell.h"

#include "gtest/gtest.h"

#ifdef __MPI
#include <mpi.h>
#endif

#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
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
        ucell.iat2it.resize(ucell.nat);
        ucell.iat2ia.resize(ucell.nat);
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
                       0,                      // init_occ_mat
                       0,                      // mixing_dftu
                       DFTU_BASE::SOC_LAYOUT_PAULI);
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

/// append_ion_step_snapshot must record an "N/A" placeholder instead of a
/// silent zero matrix when the occupation matrix does not exist yet (the PW
/// istep 0 / iter 1 case), and the real matrix body when it does.
TEST_F(DFTUBaseTest, AppendSnapshotNAPlaceholderAndReady)
{
    Plus_U_Base dftu;
    this->init_dftu(dftu, false);

    const DFTU_BASE::OccmatOutputCfg cfg = {1, 1, 5, true, 1};
    const std::string out_dir = "./";

    // fresh run: no occupation matrix has been computed or loaded
    ASSERT_FALSE(dftu.is_occmat_ready());
    DFTU_BASE::append_ion_step_snapshot(dftu,
                                        ucell,
                                        out_dir,
                                        2, // nspin
                                        1, // npol
                                        0, // istep -> occ_matg1.txt
                                        1, // iter
                                        false,
                                        false, // occmat_ready
                                        1e-6,
                                        0.5,
                                        cfg,
                                        DFTU_BASE::SOC_LAYOUT_PAULI);

    std::ifstream ifs("./occ_matg1.txt");
    ASSERT_TRUE(ifs.is_open());
    std::stringstream ss;
    ss << ifs.rdbuf();
    ifs.close();
    const std::string content = ss.str();
    EXPECT_NE(content.find("# Electronic step 1"), std::string::npos);
    EXPECT_NE(content.find("# scf_thr 1.00000000e-06"), std::string::npos);
    EXPECT_NE(content.find("# drho 5.00000000e-01"), std::string::npos);
    EXPECT_NE(content.find("\n N/A\n"), std::string::npos);
    EXPECT_EQ(content.find("Atom"), std::string::npos);

    // ready case: the matrix body is written for the atom (fresh file g2)
    DFTU_BASE::append_ion_step_snapshot(dftu,
                                        ucell,
                                        out_dir,
                                        2, // nspin
                                        1, // npol
                                        1, // istep -> occ_matg2.txt
                                        2, // iter
                                        false,
                                        true, // occmat_ready
                                        1e-6,
                                        0.5,
                                        cfg,
                                        DFTU_BASE::SOC_LAYOUT_PAULI);

    std::ifstream ifs2("./occ_matg2.txt");
    ASSERT_TRUE(ifs2.is_open());
    std::stringstream ss2;
    ss2 << ifs2.rdbuf();
    ifs2.close();
    const std::string content2 = ss2.str();
    EXPECT_NE(content2.find("Fe Atom 1 L 2"), std::string::npos);
    EXPECT_EQ(content2.find("N/A"), std::string::npos);

    std::remove("./occ_matg1.txt");
    std::remove("./occ_matg2.txt");
}

/// out_occ_mat = false must suppress both the numbered snapshot file and the
/// latest occ_mat.txt, even when the frequency gates would trigger.
TEST_F(DFTUBaseTest, OccMatSwitchDisabledWritesNothing)
{
    Plus_U_Base dftu;
    this->init_dftu(dftu, false);

    const DFTU_BASE::OccmatOutputCfg cfg = {1, 1, 5, false, 1};
    const std::string out_dir = "./";

    DFTU_BASE::append_ion_step_snapshot(dftu,
                                        ucell,
                                        out_dir,
                                        2, // nspin
                                        1, // npol
                                        0, // istep (an output ionic step)
                                        1, // iter
                                        false,
                                        true, // occmat_ready
                                        1e-6,
                                        0.5,
                                        cfg,
                                        DFTU_BASE::SOC_LAYOUT_PAULI);

    std::ifstream ifs("./occ_matg1.txt");
    EXPECT_FALSE(ifs.is_open());

    DFTU_BASE::write_latest_occmat(dftu,
                                   ucell,
                                   out_dir,
                                   2, // nspin
                                   1, // npol
                                   0, // istep
                                   2, // iter
                                   1e-6,
                                   0.5,
                                   cfg,
                                   DFTU_BASE::SOC_LAYOUT_PAULI);

    std::ifstream ifs_latest("./occ_mat.txt");
    EXPECT_FALSE(ifs_latest.is_open());
}

/// dft_plus_u = 0 (no DFT+U) must suppress both the numbered snapshot file
/// and the latest occ_mat.txt, even when out_occ_mat = true (the default).
/// This is the regression test for the abacuslite segfault: a default-
/// constructed Plus_U_Base has an empty l_channel vector, so without this
/// guard write_occup_m() dereferences a null l_channel.data() inside
/// has_l_channel(). The guard also honours the documented contract that
/// out_occ_mat only takes effect for DFT+U calculations (dft_plus_u > 0).
TEST_F(DFTUBaseTest, DftPlusUDisabledWritesNothing)
{
    Plus_U_Base dftu;  // default-constructed: l_channel is empty, mirroring
                       // the LCAO/PW esolver path when dft_plus_u == 0
    // Deliberately skip init_dftu(): the bug is that the IO functions must
    // not even reach has_l_channel() when dft_plus_u == 0.

    const DFTU_BASE::OccmatOutputCfg cfg = {1, 1, 5, true, 0};
    const std::string out_dir = "./";

    // Remove stale files so the existence check is meaningful.
    std::remove("./occ_matg1.txt");
    std::remove("./occ_mat.txt");

    DFTU_BASE::append_ion_step_snapshot(dftu,
                                        ucell,
                                        out_dir,
                                        2, // nspin
                                        1, // npol
                                        0, // istep (an output ionic step)
                                        1, // iter
                                        false,
                                        true, // occmat_ready
                                        1e-6,
                                        0.5,
                                        cfg,
                                        DFTU_BASE::SOC_LAYOUT_PAULI);

    std::ifstream ifs("./occ_matg1.txt");
    EXPECT_FALSE(ifs.is_open());

    DFTU_BASE::write_latest_occmat(dftu,
                                   ucell,
                                   out_dir,
                                   2, // nspin
                                   1, // npol
                                   0, // istep
                                   2, // iter
                                   1e-6,
                                   0.5,
                                   cfg,
                                   DFTU_BASE::SOC_LAYOUT_PAULI);

    std::ifstream ifs_latest("./occ_mat.txt");
    EXPECT_FALSE(ifs_latest.is_open());
}

/// find_first_existing_file must return the first candidate that exists,
/// in declaration order, and an empty string when none exist.
TEST(FindFirstExistingFileTest, ReturnsFirstExistingCandidate)
{
    // Use the gtest-managed temp dir for writable fixtures, and a
    // non-existent subdir for the "no candidate" case. The test never
    // creates or deletes directories itself (AGENTS.md rule 17).
    const std::string td = testing::TempDir();

    // No candidates exist -> empty string. All three names are unique to
    // this test, so they should not be present in td.
    {
        const std::vector<std::string> candidates = {"ffe_absent_1.txt",
                                                      "ffe_absent_2.txt",
                                                      "ffe_absent_3.txt"};
        EXPECT_TRUE(DFTU_BASE::find_first_existing_file(td, candidates).empty());
    }

    // Only the second candidate exists -> it is returned even though the
    // first comes earlier in the list.
    {
        const std::string fn = td + "ffe_second_only.txt";
        std::ofstream(fn).close();
        const std::vector<std::string> candidates = {"ffe_absent_1.txt",
                                                      "ffe_second_only.txt",
                                                      "ffe_absent_3.txt"};
        EXPECT_EQ(DFTU_BASE::find_first_existing_file(td, candidates), fn);
    }

    // The first candidate exists -> it takes precedence over the second.
    {
        const std::string fn = td + "ffe_first_wins.txt";
        std::ofstream(fn).close();
        const std::vector<std::string> candidates = {"ffe_first_wins.txt",
                                                      "ffe_second_only.txt",
                                                      "ffe_absent_3.txt"};
        EXPECT_EQ(DFTU_BASE::find_first_existing_file(td, candidates), fn);
    }
}

/// For init_occ_mat=2 the occupation-matrix file is read exactly once.
/// Later init_base() calls (one per ionic step in a relax run) must keep
/// the in-memory matrix instead of reading the file again, so pointing
/// the second call at a non-existent readin dir must not abort the run.
TEST_F(DFTUBaseTest, InitBaseReadsOccMatFileOnlyOnce)
{
    // Use the gtest-managed temporary directory so the test never creates
    // or deletes directories itself (see AGENTS.md rule 17).
    const std::string dir = testing::TempDir();
    const std::string fn = dir + "occ_mat.txt";

    // One Fe atom, L=2, two spin channels of 5x5, filled with distinct
    // constants so a re-read would be easy to distinguish from a
    // preserved in-memory matrix.
    {
        std::ofstream ofs(fn);
        ASSERT_TRUE(ofs.is_open());
        ofs << "# compact test fixture\n";
        ofs << " Fe Atom 1 L 2 mag 0.0\n";
        ofs << " spin 1 nelec 0.5\n";
        for (int m0 = 0; m0 < 5; ++m0)
        {
            for (int m1 = 0; m1 < 5; ++m1)
            {
                ofs << " 0.1";
            }
            ofs << "\n";
        }
        ofs << " spin 2 nelec 0.5\n";
        for (int m0 = 0; m0 < 5; ++m0)
        {
            for (int m1 = 0; m1 < 5; ++m1)
            {
                ofs << " 0.2";
            }
            ofs << "\n";
        }
    }

    Plus_U_Base dftu;
    const std::vector<int> l_channel = {2};
    const std::vector<double> hubbard_u = {0.0};
    auto call_init = [&](const std::string& readin_dir)
    {
        dftu.init_base(ucell,
                       1,                // npol
                       2,                // nspin
                       l_channel,
                       false,            // yukawa_potential
                       0.5,              // yukawa_lambda
                       readin_dir,       // global_readin_dir
                       "",               // global_out_dir
                       "none",           // init_chg
                       "cpu",            // device
                       hubbard_u,
                       0.0,              // uramping
                       2,                // init_occ_mat
                       0,                // mixing_dftu
                       DFTU_BASE::SOC_LAYOUT_PAULI);
    };

    // First ionic step: the file is read.
    call_init(dir);
    ASSERT_TRUE(dftu.is_occmat_ready());
    EXPECT_NEAR(dftu.occmat().get(0, 2, 0, 0, 0), 0.1, 1e-12);
    EXPECT_NEAR(dftu.occmat().get(0, 2, 1, 4, 4), 0.2, 1e-12);

    // Simulate a later ionic step: the in-memory matrix is preserved and
    // the file is not consulted again, so a non-existent readin dir is fine.
    call_init(dir + "does_not_exist/");
    EXPECT_TRUE(dftu.is_occmat_ready());
    EXPECT_NEAR(dftu.occmat().get(0, 2, 0, 0, 0), 0.1, 1e-12);
    EXPECT_NEAR(dftu.occmat().get(0, 2, 1, 4, 4), 0.2, 1e-12);
}

/// Reading an nspin=4 occupation-matrix file with SOC_LAYOUT_PAULI must
/// reconstruct the 4 contiguous Pauli blocks [b0, b1, b2, b3] in the 2m x 2m
/// flat buffer. The bug fixed in commit 1787365b3 was that the Im(n_ud)
/// block ("spin 12 im") was discarded, so b2 came out zero/garbage and the
/// subsequent SCF diverged. This test writes a fixture with a non-zero Im
/// block and verifies all 4 blocks at the correct flat-buffer offsets.
TEST_F(DFTUBaseTest, InitBaseReadsOccMatSocPauliRoundtrip)
{
    // Use the gtest-managed temporary directory (AGENTS.md rule 17).
    const std::string dir = testing::TempDir();
    const std::string fn = dir + "occ_mat.txt";

    // One Fe atom, L=2 -> nm=5. The 2m x 2m flat buffer (100 elements)
    // holds 4 contiguous m^2 blocks: [b0, b1, b2, b3]. The file stores
    // (n_uu, Re(n_ud), Im(n_ud), n_dd); the reader must reconstruct
    //   b0 = n_uu + n_dd
    //   b1 = 2 * Re(n_ud)
    //   b2 = 2 * Im(n_ud)   <- discarded by the buggy reader
    //   b3 = n_uu - n_dd
    const int nm = 5;
    const int m2 = nm * nm;
    std::vector<double> uu(m2), re(m2), im(m2), dd(m2);
    for (int k = 0; k < m2; ++k)
    {
        uu[k] = 1.0 + 0.01 * k;
        re[k] = 0.1 + 0.01 * k;
        im[k] = 0.01 + 0.001 * k;  // non-zero -- the bug discarded this block
        dd[k] = 0.5 + 0.01 * k;
    }

    {
        std::ofstream ofs(fn);
        ASSERT_TRUE(ofs.is_open());
        ofs << " Fe Atom 1 L 2 mag 0.0 0.0 0.0\n";
        ofs << " spin 1 nelec 0.5\n";
        for (int m0 = 0; m0 < nm; ++m0)
        {
            for (int m1 = 0; m1 < nm; ++m1)
                ofs << " " << uu[m0 * nm + m1];
            ofs << "\n";
        }
        ofs << " spin 12 re\n";
        for (int m0 = 0; m0 < nm; ++m0)
        {
            for (int m1 = 0; m1 < nm; ++m1)
                ofs << " " << re[m0 * nm + m1];
            ofs << "\n";
        }
        ofs << " spin 12 im\n";
        for (int m0 = 0; m0 < nm; ++m0)
        {
            for (int m1 = 0; m1 < nm; ++m1)
                ofs << " " << im[m0 * nm + m1];
            ofs << "\n";
        }
        ofs << " spin 2 nelec 0.5\n";
        for (int m0 = 0; m0 < nm; ++m0)
        {
            for (int m1 = 0; m1 < nm; ++m1)
                ofs << " " << dd[m0 * nm + m1];
            ofs << "\n";
        }
    }

    Plus_U_Base dftu;
    const std::vector<int> l_channel = {2};
    const std::vector<double> hubbard_u = {0.0};
    dftu.init_base(ucell,
                   2,                // npol (nspin=4 requires npol=2)
                   4,                // nspin
                   l_channel,
                   false,            // yukawa_potential
                   0.5,              // yukawa_lambda
                   dir,              // global_readin_dir
                   "",               // global_out_dir
                   "none",           // init_chg
                   "cpu",            // device
                   hubbard_u,
                   0.0,              // uramping
                   2,                // init_occ_mat
                   0,                // mixing_dftu
                   DFTU_BASE::SOC_LAYOUT_PAULI);

    ASSERT_TRUE(dftu.is_occmat_ready());

    // The 2m x 2m flat buffer holds 4 contiguous m^2 blocks: [b0, b1, b2, b3].
    const ModuleBase::matrix& occ0 = dftu.occmat().mat(0, 2, 0);
    ASSERT_EQ(occ0.nr * occ0.nc, 4 * m2);

    for (int k = 0; k < m2; ++k)
    {
        EXPECT_NEAR(occ0.c[0 * m2 + k], uu[k] + dd[k], 1e-12)    // b0
            << "Pauli block b0 mismatch at k=" << k;
        EXPECT_NEAR(occ0.c[1 * m2 + k], 2.0 * re[k], 1e-12)      // b1
            << "Pauli block b1 mismatch at k=" << k;
        EXPECT_NEAR(occ0.c[2 * m2 + k], 2.0 * im[k], 1e-12)      // b2 -- the bug
            << "Pauli block b2 mismatch at k=" << k;
        EXPECT_NEAR(occ0.c[3 * m2 + k], uu[k] - dd[k], 1e-12)    // b3
            << "Pauli block b3 mismatch at k=" << k;
    }
}

// Reading the occupation-matrix file broadcasts on MPI_COMM_WORLD, so
// the test binary must initialize MPI even when ctest launches it as a
// single process.
int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
#endif
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    MPI_Finalize();
#endif
    return result;
}
