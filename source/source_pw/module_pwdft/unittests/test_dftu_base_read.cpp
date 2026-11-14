#include "dftu_base_test_fixture.h"
#include "source_pw/module_pwdft/dftu_base.h"
#include "source_pw/module_pwdft/dftu_pw_tools.h"
#include "source_cell/atom_spec.h"

#ifdef __MPI
#include <mpi.h>
#endif

#include <cstdio>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <numeric>


namespace
{

void write_constant_matrix(std::ostream& os,
                           const int size,
                           const double value,
                           const int first_row_columns,
                           const bool inline_comment = false)
{
    for (int row = 0; row < size; ++row)
    {
        const int columns = (row == 0 && first_row_columns >= 0)
                                ? first_row_columns
                                : size;
        for (int column = 0; column < columns; ++column)
        {
            os << " " << value;
        }
        if (inline_comment)
        {
            os << " # matrix row";
        }
        os << "\n";
    }
}

void write_compact_collinear_body(std::ostream& os)
{
    os << "Fe Atom 1 L 2 mag 0.0\n";
    os << "spin 1 nelec 0.5\n";
    write_constant_matrix(os, 5, 0.1, 5);
    os << "spin 2 nelec 0.5\n";
    write_constant_matrix(os, 5, 0.2, 5);
}

void write_compact_collinear_atom(std::ostream& os,
                                  const int atom_index,
                                  const double spin_up,
                                  const double spin_down)
{
    os << "Fe Atom " << atom_index << " L 2 mag 0.0\n";
    os << "spin 1 nelec 0.5\n";
    write_constant_matrix(os, 5, spin_up, 5);
    os << "spin 2 nelec 0.5\n";
    write_constant_matrix(os, 5, spin_down, 5);
}

} // namespace


/// init_occ_mat=2 loads the fixed occupation matrix during run initialization.
TEST_F(DFTUBaseTest, InitBaseReadsFixedOccMat)
{
    // Use the gtest-managed temporary directory so the test never creates
    // or deletes directories itself (see AGENTS.md rule 17).
    const std::string dir = testing::TempDir();
    const std::string fn = dir + "occ_mat.txt";

    // One Fe atom, L=2, two spin channels of 5x5, filled with distinct
    // constants.
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
    dftu.init_base(ucell,
                   1,                // npol
                   2,                // nspin
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

TEST_F(DFTUBaseTest, ReadOccMatCompactHandlesLongComment)
{
    const std::string fn = testing::TempDir() + "compact_long_comment.txt";
    {
        std::ofstream ofs(fn);
        ASSERT_TRUE(ofs.is_open());
        ofs << "# " << std::string(300, 'x') << "\n";
        ofs << "Fe Atom 1 L 2 mag 0.0 # atom header\n";
        ofs << "spin 1 nelec 0.5 # spin header\n";
        write_constant_matrix(ofs, 5, 0.1, 5, true);
        ofs << "spin 2 nelec 0.5 # spin header\n";
        write_constant_matrix(ofs, 5, 0.2, 5, true);
    }

    const std::vector<int> l_channel = {2};
    OccupationMatrix occ;
    occ.init(ucell, l_channel, 2, 1);
    DFTU_BASE::read_occup_m(ucell,
                            occ,
                            l_channel,
                            fn,
                            2,
                            1,
                            DFTU_BASE::SOC_LAYOUT_PAULI);
    EXPECT_DOUBLE_EQ(occ.get(0, 2, 0, 0, 0), 0.1);
    EXPECT_DOUBLE_EQ(occ.get(0, 2, 1, 4, 4), 0.2);
}

TEST_F(DFTUBaseTest, ReadOccMatCompactMultipleAtomsSkipsUncorrelatedSpecies)
{
    UnitCell multi_cell;
    multi_cell.ntype = 2;
    multi_cell.nat = 4;
    multi_cell.atoms = new Atom[multi_cell.ntype];
    multi_cell.set_atom_flag = true;

    multi_cell.atoms[0].label = "O";
    multi_cell.atoms[0].na = 2;
    multi_cell.atoms[0].nwl = 1;
    multi_cell.atoms[0].nw = 4;

    Atom& fe = multi_cell.atoms[1];
    fe.label = "Fe";
    fe.na = 2;
    fe.nwl = 2;
    fe.nw = 9;
    fe.iw2l.resize(fe.nw);
    fe.iw2n.resize(fe.nw);
    fe.iw2m.resize(fe.nw);
    int iw = 0;
    for (int l = 0; l <= fe.nwl; ++l)
    {
        for (int m = 0; m < 2 * l + 1; ++m)
        {
            fe.iw2l[iw] = l;
            fe.iw2n[iw] = 0;
            fe.iw2m[iw] = m;
            ++iw;
        }
    }

    multi_cell.set_iat2itia();
    multi_cell.itia2iat.create(multi_cell.ntype, 2);
    multi_cell.itia2iat(0, 0) = 0;
    multi_cell.itia2iat(0, 1) = 1;
    multi_cell.itia2iat(1, 0) = 2;
    multi_cell.itia2iat(1, 1) = 3;
    multi_cell.set_iat2iwt(1);

    const std::string fn = testing::TempDir() + "compact_multi_atom.txt";
    {
        std::ofstream ofs(fn);
        ASSERT_TRUE(ofs.is_open());
        write_compact_collinear_atom(ofs, 3, 0.1, 0.2);
        write_compact_collinear_atom(ofs, 4, 0.3, 0.4);
    }

    const std::vector<int> l_channel = {-1, 2};
    OccupationMatrix occ;
    occ.init(multi_cell, l_channel, 2, 1);
    DFTU_BASE::read_occup_m(multi_cell,
                            occ,
                            l_channel,
                            fn,
                            2,
                            1,
                            DFTU_BASE::SOC_LAYOUT_PAULI);

    EXPECT_DOUBLE_EQ(occ.get(2, 2, 0, 0, 0), 0.1);
    EXPECT_DOUBLE_EQ(occ.get(2, 2, 1, 4, 4), 0.2);
    EXPECT_DOUBLE_EQ(occ.get(3, 2, 0, 0, 0), 0.3);
    EXPECT_DOUBLE_EQ(occ.get(3, 2, 1, 4, 4), 0.4);
}

TEST_F(DFTUBaseTest, ReadOccMatCompactSocRealSpinBasis)
{
    const std::string fn = testing::TempDir() + "compact_soc_real.txt";
    {
        std::ofstream ofs(fn);
        ASSERT_TRUE(ofs.is_open());
        ofs << "Fe Atom 1 L 2 mag 0.0 0.0 0.0\n";
        ofs << "spin 1 nelec 0.5\n";
        write_constant_matrix(ofs, 5, 1.0, 5);
        ofs << "spin 12 re\n";
        write_constant_matrix(ofs, 5, 2.0, 5);
        ofs << "spin 12 im\n";
        write_constant_matrix(ofs, 5, 3.0, 5);
        ofs << "spin 2 nelec 0.5\n";
        write_constant_matrix(ofs, 5, 4.0, 5);
    }

    const std::vector<int> l_channel = {2};
    OccupationMatrix occ;
    occ.init(ucell, l_channel, 4, 2);
    DFTU_BASE::read_occup_m(ucell,
                            occ,
                            l_channel,
                            fn,
                            4,
                            2,
                            DFTU_BASE::SOC_LAYOUT_SPIN_BASIS_REAL);
    EXPECT_DOUBLE_EQ(occ.get(0, 2, 0, 0, 0), 1.0);
    EXPECT_DOUBLE_EQ(occ.get(0, 2, 0, 0, 5), 2.0);
    EXPECT_DOUBLE_EQ(occ.get(0, 2, 0, 5, 0), 2.0);
    EXPECT_DOUBLE_EQ(occ.get(0, 2, 0, 5, 5), 4.0);
}

TEST_F(DFTUBaseTest, ReadOccMatRejectsLegacyHeader)
{
    const std::string fn = testing::TempDir() + "legacy_header.txt";
    {
        std::ofstream ofs(fn);
        ASSERT_TRUE(ofs.is_open());
        ofs << "Fe 1 Atom=1 L=2 0 0 0\n";
    }
    const std::vector<int> l_channel = {2};
    OccupationMatrix occ;
    occ.init(ucell, l_channel, 2, 1);
    EXPECT_EXIT(DFTU_BASE::read_occup_m(ucell,
                                       occ,
                                       l_channel,
                                       fn,
                                       2,
                                       1,
                                       DFTU_BASE::SOC_LAYOUT_PAULI),
                testing::ExitedWithCode(1),
                "");
}

TEST_F(DFTUBaseTest, ReadOccMatRejectsShortMatrixRow)
{
    const std::string fn = testing::TempDir() + "short_matrix_row.txt";
    {
        std::ofstream ofs(fn);
        ASSERT_TRUE(ofs.is_open());
        ofs << "Fe Atom 1 L 2 mag 0.0\n";
        ofs << "spin 1 nelec 0.5\n";
        write_constant_matrix(ofs, 5, 0.1, 4);
        ofs << "spin 2 nelec 0.5\n";
        write_constant_matrix(ofs, 5, 0.2, 5);
    }
    const std::vector<int> l_channel = {2};
    OccupationMatrix occ;
    occ.init(ucell, l_channel, 2, 1);
    EXPECT_EXIT(DFTU_BASE::read_occup_m(ucell,
                                       occ,
                                       l_channel,
                                       fn,
                                       2,
                                       1,
                                       DFTU_BASE::SOC_LAYOUT_PAULI),
                testing::ExitedWithCode(1),
                "");
}

TEST_F(DFTUBaseTest, ReadOccMatRejectsLongMatrixRow)
{
    const std::string fn = testing::TempDir() + "long_matrix_row.txt";
    {
        std::ofstream ofs(fn);
        ASSERT_TRUE(ofs.is_open());
        ofs << "Fe Atom 1 L 2 mag 0.0\n";
        ofs << "spin 1 nelec 0.5\n";
        write_constant_matrix(ofs, 5, 0.1, 6);
        ofs << "spin 2 nelec 0.5\n";
        write_constant_matrix(ofs, 5, 0.2, 5);
    }
    const std::vector<int> l_channel = {2};
    OccupationMatrix occ;
    occ.init(ucell, l_channel, 2, 1);
    EXPECT_EXIT(DFTU_BASE::read_occup_m(ucell,
                                       occ,
                                       l_channel,
                                       fn,
                                       2,
                                       1,
                                       DFTU_BASE::SOC_LAYOUT_PAULI),
                testing::ExitedWithCode(1),
                "");
}

TEST_F(DFTUBaseTest, ReadOccMatRejectsMissingSocBlock)
{
    const std::string fn = testing::TempDir() + "missing_soc_block.txt";
    {
        std::ofstream ofs(fn);
        ASSERT_TRUE(ofs.is_open());
        ofs << "Fe Atom 1 L 2 mag 0.0 0.0 0.0\n";
        ofs << "spin 1 nelec 0.5\n";
        write_constant_matrix(ofs, 5, 1.0, 5);
        ofs << "spin 12 re\n";
        write_constant_matrix(ofs, 5, 2.0, 5);
        ofs << "spin 2 nelec 0.5\n";
    }
    const std::vector<int> l_channel = {2};
    OccupationMatrix occ;
    occ.init(ucell, l_channel, 4, 2);
    EXPECT_EXIT(DFTU_BASE::read_occup_m(ucell,
                                       occ,
                                       l_channel,
                                       fn,
                                       4,
                                       2,
                                       DFTU_BASE::SOC_LAYOUT_PAULI),
                testing::ExitedWithCode(1),
                "");
}

TEST_F(DFTUBaseTest, ReadOccMatRejectsExtraContent)
{
    const std::string fn = testing::TempDir() + "extra_content.txt";
    {
        std::ofstream ofs(fn);
        ASSERT_TRUE(ofs.is_open());
        write_compact_collinear_body(ofs);
        ofs << "unexpected_trailing_token_that_is_longer_than_twenty_characters\n";
    }
    const std::vector<int> l_channel = {2};
    OccupationMatrix occ;
    occ.init(ucell, l_channel, 2, 1);
    EXPECT_EXIT(DFTU_BASE::read_occup_m(ucell,
                                       occ,
                                       l_channel,
                                       fn,
                                       2,
                                       1,
                                       DFTU_BASE::SOC_LAYOUT_PAULI),
                testing::ExitedWithCode(1),
                "");
}
