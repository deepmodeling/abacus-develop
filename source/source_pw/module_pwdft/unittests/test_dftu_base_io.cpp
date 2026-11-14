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

TEST_F(DFTUBaseTest, CompactWriterReaderCollinearRoundtrip)
{
    Plus_U_Base source;
    init_dftu_for_spin(source, false, 2, 1, DFTU_BASE::SOC_LAYOUT_PAULI);
    for (int spin = 0; spin < 2; ++spin)
    {
        for (int row = 0; row < 5; ++row)
        {
            for (int column = 0; column < 5; ++column)
            {
                const double value = 0.1 * (spin + 1) + 0.01 * row + 0.001 * column;
                source.occmat().set(0, 2, spin, row, column, value);
            }
        }
    }

    const std::string filename = testing::TempDir() + "compact_writer_collinear.txt";
    {
        std::ofstream output(filename.c_str());
        ASSERT_TRUE(output.is_open());
        DFTU_BASE::write_occup_m(source,
                                 ucell,
                                 output,
                                 DFTU_BASE::SOC_LAYOUT_PAULI);
    }

    OccupationMatrix restored;
    const std::vector<int> l_channel = {2};
    restored.init(ucell, l_channel, 2, 1);
    DFTU_BASE::read_occup_m(ucell,
                            restored,
                            l_channel,
                            filename,
                            2,
                            1,
                            DFTU_BASE::SOC_LAYOUT_PAULI);

    for (int spin = 0; spin < 2; ++spin)
    {
        for (int row = 0; row < 5; ++row)
        {
            for (int column = 0; column < 5; ++column)
            {
                EXPECT_NEAR(restored.get(0, 2, spin, row, column),
                            source.occmat().get(0, 2, spin, row, column),
                            1e-8);
            }
        }
    }

    // Assert the diagnostic header values that read_occup_m only validates
    // for finiteness: per-spin nelec (trace) and collinear mag (up - down).
    double expected_nelec[2] = {0.0, 0.0};
    for (int spin = 0; spin < 2; ++spin)
    {
        for (int m = 0; m < 5; ++m)
        {
            expected_nelec[spin] += source.occmat().get(0, 2, spin, m, m);
        }
    }
    const double expected_mag = expected_nelec[0] - expected_nelec[1];

    std::ifstream input(filename.c_str());
    ASSERT_TRUE(input.is_open());
    std::string token;
    double parsed_mag = 0.0;
    double parsed_nelec[2] = {0.0, 0.0};
    int spin_index = -1;
    while (input >> token)
    {
        if (token == "mag")
        {
            ASSERT_TRUE(input >> parsed_mag);
        }
        else if (token == "spin")
        {
            ASSERT_TRUE(input >> spin_index);
            ASSERT_GE(spin_index, 1);
            ASSERT_LE(spin_index, 2);
        }
        else if (token == "nelec")
        {
            ASSERT_EQ(spin_index == 1 || spin_index == 2, true);
            ASSERT_TRUE(input >> parsed_nelec[spin_index - 1]);
        }
    }
    EXPECT_NEAR(parsed_mag, expected_mag, 1e-7);
    EXPECT_NEAR(parsed_nelec[0], expected_nelec[0], 1e-7);
    EXPECT_NEAR(parsed_nelec[1], expected_nelec[1], 1e-7);
}

TEST_F(DFTUBaseTest, CompactWriterReaderSocPauliRoundtrip)
{
    Plus_U_Base source;
    const DFTU_BASE::OccmatSocLayout layout = DFTU_BASE::lcao_soc_layout(1);
    EXPECT_EQ(layout, DFTU_BASE::SOC_LAYOUT_PAULI);
    init_dftu_for_spin(source, false, 4, 2, layout);
    const int block_size = 25;
    for (int block = 0; block < 4; ++block)
    {
        for (int index = 0; index < block_size; ++index)
        {
            source.occmat().mat(0, 2, 0).c[block * block_size + index]
                = 0.1 * (block + 1) + 0.001 * index;
        }
    }

    const std::string filename = testing::TempDir() + "compact_writer_soc_pauli.txt";
    {
        std::ofstream output(filename.c_str());
        ASSERT_TRUE(output.is_open());
        DFTU_BASE::write_occup_m(source, ucell, output, layout);
    }

    OccupationMatrix restored;
    const std::vector<int> l_channel = {2};
    restored.init(ucell, l_channel, 4, 2);
    DFTU_BASE::read_occup_m(ucell,
                            restored,
                            l_channel,
                            filename,
                            4,
                            2,
                            layout);

    for (int index = 0; index < 4 * block_size; ++index)
    {
        EXPECT_NEAR(restored.mat(0, 2, 0).c[index],
                    source.occmat().mat(0, 2, 0).c[index],
                    1e-8);
    }

    // Assert the SOC diagnostic header: mag x/y/z and per-spin nelec.
    // Pauli blocks are charge/spin_x/spin_y/spin_z, so the writer reports
    // mag_x = sum(spin_x diag), mag_y = sum(spin_y diag), mag_z = sum(spin_z diag)
    // and nelec = 0.5 * sum(charge +/- spin_z) on the diagonal.
    const ModuleBase::matrix& occ = source.occmat().mat(0, 2, 0);
    double expected_mag[3] = {0.0, 0.0, 0.0};
    double expected_nelec_soc[2] = {0.0, 0.0};
    for (int m = 0; m < 5; ++m)
    {
        const int diag = m * 5 + m;
        expected_mag[0] += occ.c[1 * block_size + diag];
        expected_mag[1] += occ.c[2 * block_size + diag];
        expected_mag[2] += occ.c[3 * block_size + diag];
        expected_nelec_soc[0] += 0.5 * (occ.c[0 * block_size + diag] + occ.c[3 * block_size + diag]);
        expected_nelec_soc[1] += 0.5 * (occ.c[0 * block_size + diag] - occ.c[3 * block_size + diag]);
    }
    // The fixture fills every block entry, including a nonzero spin_y block,
    // so mag_y must be nonzero for this assertion to be meaningful.
    EXPECT_NE(expected_mag[1], 0.0);

    std::ifstream soc_input(filename.c_str());
    ASSERT_TRUE(soc_input.is_open());
    std::string token;
    double parsed_mag_soc[3] = {0.0, 0.0, 0.0};
    double parsed_nelec_soc[2] = {0.0, 0.0};
    int spin_index = -1;
    while (soc_input >> token)
    {
        if (token == "mag")
        {
            ASSERT_TRUE(soc_input >> parsed_mag_soc[0] >> parsed_mag_soc[1] >> parsed_mag_soc[2]);
        }
        else if (token == "spin")
        {
            ASSERT_TRUE(soc_input >> spin_index);
        }
        else if (token == "nelec")
        {
            ASSERT_EQ(spin_index == 1 || spin_index == 2, true);
            ASSERT_TRUE(soc_input >> parsed_nelec_soc[spin_index - 1]);
        }
    }
    EXPECT_NEAR(parsed_mag_soc[0], expected_mag[0], 1e-7);
    EXPECT_NEAR(parsed_mag_soc[1], expected_mag[1], 1e-7);
    EXPECT_NEAR(parsed_mag_soc[2], expected_mag[2], 1e-7);
    EXPECT_NEAR(parsed_nelec_soc[0], expected_nelec_soc[0], 1e-7);
    EXPECT_NEAR(parsed_nelec_soc[1], expected_nelec_soc[1], 1e-7);
}

TEST_F(DFTUBaseTest, CompactWriterReaderSocRealSpinBasisRoundtrip)
{
    Plus_U_Base source;
    const DFTU_BASE::OccmatSocLayout layout = DFTU_BASE::lcao_soc_layout(2);
    EXPECT_EQ(layout, DFTU_BASE::SOC_LAYOUT_SPIN_BASIS_REAL);
    init_dftu_for_spin(source, false, 4, 2, layout);
    const int orbital_size = 5;
    for (int row = 0; row < orbital_size; ++row)
    {
        for (int column = 0; column < orbital_size; ++column)
        {
            const double up_up = 0.1 + 0.01 * row + 0.001 * column;
            const double up_down = 0.2 + 0.01 * row + 0.001 * column;
            const double down_up = 0.2 + 0.01 * column + 0.001 * row;
            const double down_down = 0.3 + 0.01 * row + 0.001 * column;
            source.occmat().set(0, 2, 0, row, column, up_up);
            source.occmat().set(0, 2, 0, row, column + orbital_size, up_down);
            source.occmat().set(0, 2, 0, row + orbital_size, column, down_up);
            source.occmat().set(0, 2, 0, row + orbital_size, column + orbital_size, down_down);
        }
    }

    const std::string filename = testing::TempDir() + "compact_writer_soc_real.txt";
    {
        std::ofstream output(filename.c_str());
        ASSERT_TRUE(output.is_open());
        DFTU_BASE::write_occup_m(source, ucell, output, layout);
    }

    OccupationMatrix restored;
    const std::vector<int> l_channel = {2};
    restored.init(ucell, l_channel, 4, 2);
    DFTU_BASE::read_occup_m(ucell,
                            restored,
                            l_channel,
                            filename,
                            4,
                            2,
                            layout);

    for (int row = 0; row < 2 * orbital_size; ++row)
    {
        for (int column = 0; column < 2 * orbital_size; ++column)
        {
            EXPECT_NEAR(restored.get(0, 2, 0, row, column),
                        source.occmat().get(0, 2, 0, row, column),
                        1e-8);
        }
    }
}
