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


/// append_ion_step_snapshot must record an "N/A" placeholder instead of a
/// silent zero matrix when the occupation matrix does not exist yet (the PW
/// istep 0 / iter 1 case), and the real matrix body when it does.
TEST_F(DFTUBaseTest, AppendSnapshotNAPlaceholderAndReady)
{
    Plus_U_Base dftu;
    this->init_dftu(dftu, false);

    const DFTU_BASE::OccmatOutputCfg cfg = {1, 1, 5, true, 1};
    const std::string out_dir = testing::TempDir();
    const DFTU_BASE::OccmatStepInfo first_step = {0, 1, 1e-6, 0.5};
    const DFTU_BASE::OccmatStepInfo second_step = {1, 2, 1e-6, 0.5};

    // fresh run: no occupation matrix has been computed or loaded
    ASSERT_FALSE(dftu.is_occmat_ready());
    DFTU_BASE::append_ion_step_snapshot(dftu,
                                        ucell,
                                        out_dir,
                                        first_step,
                                        false,
                                        false, // occmat_ready
                                        cfg,
                                        DFTU_BASE::SOC_LAYOUT_PAULI);

    const std::string fn1 = out_dir + "occ_matg1.txt";
    std::ifstream ifs(fn1.c_str());
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
                                        second_step,
                                        false,
                                        true, // occmat_ready
                                        cfg,
                                        DFTU_BASE::SOC_LAYOUT_PAULI);

    const std::string fn2 = out_dir + "occ_matg2.txt";
    std::ifstream ifs2(fn2.c_str());
    ASSERT_TRUE(ifs2.is_open());
    std::stringstream ss2;
    ss2 << ifs2.rdbuf();
    ifs2.close();
    const std::string content2 = ss2.str();
    EXPECT_NE(content2.find("Fe Atom 1 L 2"), std::string::npos);
    EXPECT_EQ(content2.find("N/A"), std::string::npos);
    EXPECT_EQ(content2.find("# eigenvalues"), std::string::npos);

    // Clean up so subsequent tests checking "file not written" are not
    // fooled by leftovers.  std::remove is Rule 17 compliant.
    std::remove((out_dir + "occ_matg1.txt").c_str());
    std::remove((out_dir + "occ_matg2.txt").c_str());
}

/// out_occ_mat = false must suppress both the numbered snapshot file and the
/// latest occ_mat.txt, even when the frequency gates would trigger.
TEST_F(DFTUBaseTest, OccMatSwitchDisabledWritesNothing)
{
    Plus_U_Base dftu;
    this->init_dftu(dftu, false);

    const DFTU_BASE::OccmatOutputCfg cfg = {1, 1, 5, false, 1};
    const std::string out_dir = testing::TempDir();
    const DFTU_BASE::OccmatStepInfo step = {0, 1, 1e-6, 0.5};

    // Remove stale files from previous test binary runs so the
    // "file not written" assertions are meaningful.
    std::remove((out_dir + "occ_matg1.txt").c_str());
    std::remove((out_dir + "occ_mat.txt").c_str());

    DFTU_BASE::append_ion_step_snapshot(dftu,
                                        ucell,
                                        out_dir,
                                        step,
                                        false,
                                        true, // occmat_ready
                                        cfg,
                                        DFTU_BASE::SOC_LAYOUT_PAULI);

    const std::string fn_g1 = out_dir + "occ_matg1.txt";
    std::ifstream ifs(fn_g1.c_str());
    EXPECT_FALSE(ifs.is_open());

    DFTU_BASE::write_latest_occmat(dftu,
                                   ucell,
                                   out_dir,
                                   step,
                                   cfg,
                                   DFTU_BASE::SOC_LAYOUT_PAULI);

    const std::string fn_latest = out_dir + "occ_mat.txt";
    std::ifstream ifs_latest(fn_latest.c_str());
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
    const std::string out_dir = testing::TempDir() + "guard_test/";
    const DFTU_BASE::OccmatStepInfo step = {0, 1, 1e-6, 0.5};

    DFTU_BASE::append_ion_step_snapshot(dftu,
                                        ucell,
                                        out_dir,
                                        step,
                                        false,
                                        true, // occmat_ready
                                        cfg,
                                        DFTU_BASE::SOC_LAYOUT_PAULI);

    const std::string fn_g1 = out_dir + "occ_matg1.txt";
    std::ifstream ifs(fn_g1.c_str());
    EXPECT_FALSE(ifs.is_open());

    DFTU_BASE::write_latest_occmat(dftu,
                                   ucell,
                                   out_dir,
                                   step,
                                   cfg,
                                   DFTU_BASE::SOC_LAYOUT_PAULI);

    const std::string fn_latest = out_dir + "occ_mat.txt";
    std::ifstream ifs_latest(fn_latest.c_str());
    EXPECT_FALSE(ifs_latest.is_open());
}

TEST_F(DFTUBaseTest, LatestOccMatIncludesSortedCollinearEigenvalues)
{
    Plus_U_Base dftu;
    init_dftu(dftu, false);
    for (int spin = 0; spin < 2; ++spin)
    {
        for (int m = 0; m < 5; ++m)
        {
            const double occupation = spin == 0 ? 0.8 : 0.2;
            dftu.occmat().set(0, 2, spin, m, m, occupation);
        }
    }
    dftu.occmat().set(0, 2, 0, 0, 1, 0.1);
    dftu.occmat().set(0, 2, 0, 1, 0, 0.1);

    const DFTU_BASE::OccmatOutputCfg cfg = {0, 1, 5, true, 1};
    const DFTU_BASE::OccmatStepInfo step = {0, 2, 1e-6, 1e-7};
    const std::string filename = testing::TempDir() + "occ_mat.txt";
    DFTU_BASE::write_latest_occmat(dftu, ucell, testing::TempDir(), step, cfg,
                                   DFTU_BASE::SOC_LAYOUT_PAULI);

    std::ifstream input(filename.c_str());
    ASSERT_TRUE(input.is_open());
    std::string line;
    int eigenvalue_lines = 0;
    while (std::getline(input, line))
    {
        if (line.find("# eigenvalues atom 1 L 2 spin ") != 0)
        {
            continue;
        }
        ++eigenvalue_lines;
        ASSERT_LE(eigenvalue_lines, 2);
        const std::size_t separator = line.find(':');
        ASSERT_NE(separator, std::string::npos);
        const std::size_t values_start = separator + 1;
        const std::string values_text = line.substr(values_start);
        std::istringstream values(values_text);
        const double expected[2][5] = {{0.7, 0.8, 0.8, 0.8, 0.9},
                                       {0.2, 0.2, 0.2, 0.2, 0.2}};
        double value = 0.0;
        for (int m = 0; m < 5; ++m)
        {
            ASSERT_TRUE(values >> value);
            EXPECT_NEAR(value, expected[eigenvalue_lines - 1][m], 1e-8);
        }
        EXPECT_FALSE(values >> value);
    }
    EXPECT_EQ(eigenvalue_lines, 2);

    OccupationMatrix restored;
    const std::vector<int> l_channel = {2};
    restored.init(ucell, l_channel, 2, 1);
    DFTU_BASE::read_occup_m(ucell, restored, l_channel, filename, 2, 1,
                            DFTU_BASE::SOC_LAYOUT_PAULI);
    EXPECT_NEAR(restored.get(0, 2, 0, 0, 1), 0.1, 1e-8);
    std::remove(filename.c_str());
}

TEST_F(DFTUBaseTest, SocEigenvaluesIncludeImaginarySpinCoupling)
{
    Plus_U_Base dftu;
    init_dftu_for_spin(dftu, false, 4, 2, DFTU_BASE::SOC_LAYOUT_PAULI);
    const int block_size = 25;
    for (int m = 0; m < 5; ++m)
    {
        const int index = m * 5 + m;
        dftu.occmat().mat(0, 2, 0).c[index] = 1.0;
        dftu.occmat().mat(0, 2, 0).c[2 * block_size + index] = 0.6;
        dftu.occmat().mat(0, 2, 0).c[3 * block_size + index] = 0.6;
    }
    std::ostringstream output;
    DFTU_BASE::write_occup_eigenvalues(dftu, ucell, output,
                                       DFTU_BASE::SOC_LAYOUT_PAULI);
    const std::string line = output.str();
    ASSERT_EQ(line.find("# eigenvalues atom 1 L 2 spinor:"), 0);
    const std::size_t separator = line.find(':');
    ASSERT_NE(separator, std::string::npos);
    const std::size_t values_start = separator + 1;
    const std::string values_text = line.substr(values_start);
    std::istringstream values(values_text);
    const double low = 0.5 - std::sqrt(0.18);
    const double high = 0.5 + std::sqrt(0.18);
    double value = 0.0;
    for (int i = 0; i < 10; ++i)
    {
        ASSERT_TRUE(values >> value);
        const double expected = i < 5 ? low : high;
        EXPECT_NEAR(value, expected, 1e-8);
    }
    EXPECT_FALSE(values >> value);
}

TEST_F(DFTUBaseTest, CollinearEigenvaluesStopAtJacobiCutoff)
{
    Plus_U_Base dftu;
    init_dftu(dftu, false);
    dftu.occmat().set(0, 2, 0, 0, 0, 0.5);
    dftu.occmat().set(0, 2, 0, 1, 1, 0.5);
    dftu.occmat().set(0, 2, 0, 0, 1, 1.0e-10);
    dftu.occmat().set(0, 2, 0, 1, 0, 1.0e-10);

    std::ostringstream output;
    DFTU_BASE::write_occup_eigenvalues(dftu, ucell, output,
                                       DFTU_BASE::SOC_LAYOUT_PAULI);
    EXPECT_NE(output.str().find("# eigenvalues atom 1 L 2 spin 1:"), std::string::npos);
}

TEST_F(DFTUBaseTest, NonmagneticEigenvaluesUseOnlyMaintainedChannel)
{
    Plus_U_Base dftu;
    init_dftu_for_spin(dftu, false, 1, 1, DFTU_BASE::SOC_LAYOUT_PAULI);
    dftu.occmat().set(0, 2, 0, 0, 0, 0.6);
    dftu.occmat().set(0, 2, 1, 0, 0, 0.1); // unused or stale channel

    std::ostringstream output;
    DFTU_BASE::write_occup_eigenvalues(dftu, ucell, output,
                                       DFTU_BASE::SOC_LAYOUT_PAULI);
    const std::string result = output.str();
    EXPECT_EQ(result.find("spin 2"), std::string::npos);
    EXPECT_EQ(result.find("spin 1"), std::string::npos);
    EXPECT_NE(result.find("# eigenvalues atom 1 L 2 nonmagnetic: 0.00000000 0.00000000 0.00000000 0.00000000 0.60000000"),
              std::string::npos);
    EXPECT_EQ(result.find("\n# eigenvalues", 1), std::string::npos);
}

TEST_F(DFTUBaseTest, NonmagneticPwEigenvaluesKeepSpinSummedWeight)
{
    Plus_U_Base dftu;
    init_dftu_for_spin(dftu, false, 1, 1, DFTU_BASE::SOC_LAYOUT_PAULI);
    std::vector<std::complex<double>> projector_overlap(5, {0.0, 0.0});
    projector_overlap[0] = {1.0, 0.0};
    ModuleBase::matrix band_weight(1, 1);
    band_weight(0, 0) = 2.0; // nspin=1 k-point weights include spin degeneracy
    pw::accumulate_occ_scalar(dftu.occmat().mat(0, 2, 0).c,
                              projector_overlap.data(), 1, 5, 0, 0, 5,
                              band_weight, 0);

    std::ostringstream output;
    DFTU_BASE::write_occup_eigenvalues(dftu, ucell, output,
                                       DFTU_BASE::SOC_LAYOUT_PAULI);
    EXPECT_NE(output.str().find("nonmagnetic: 0.00000000 0.00000000 0.00000000 0.00000000 2.00000000"),
              std::string::npos);
}

TEST_F(DFTUBaseTest, LcaoMethodTwoSocEigenvaluesPreserveCrossOrbitalCoupling)
{
    Plus_U_Base dftu;
    init_dftu_for_spin(dftu, false, 4, 2, DFTU_BASE::SOC_LAYOUT_SPIN_BASIS_REAL);
    const int orbital_size = 5;
    const int matrix_size = 2 * orbital_size;
    const int cross_spin_orbital = orbital_size + 1;
    for (int m = 0; m < matrix_size; ++m)
    {
        dftu.occmat().set(0, 2, 0, m, m, 0.5);
    }
    dftu.occmat().set(0, 2, 0, 0, cross_spin_orbital, 0.2);
    dftu.occmat().set(0, 2, 0, cross_spin_orbital, 0, 0.2);

    std::ostringstream output;
    DFTU_BASE::write_occup_eigenvalues(dftu, ucell, output,
                                       DFTU_BASE::SOC_LAYOUT_SPIN_BASIS_REAL);
    const std::string line = output.str();
    const std::size_t separator = line.find(':');
    ASSERT_NE(separator, std::string::npos);
    const std::size_t values_start = separator + 1;
    const std::string values_text = line.substr(values_start);
    std::istringstream values(values_text);
    double value = 0.0;
    for (int i = 0; i < matrix_size; ++i)
    {
        ASSERT_TRUE(values >> value);
        const double expected = i == 0 ? 0.3 : (i == matrix_size - 1 ? 0.7 : 0.5);
        EXPECT_NEAR(value, expected, 1e-8);
    }
    EXPECT_FALSE(values >> value);
}

/// When DFT+U is off (dft_plus_u <= 0) append_ion_step_snapshot must
/// return immediately, BEFORE the iter >= 1 validation.  This guard-first
/// ordering was introduced in step 3/3 of the EXX restart refactor so that
/// a transient iter=0 (e.g. from an upstream EXX rerun signal) does not
/// trigger WARNING_QUIT on a non-DFT+U calculation.  Passing iter=0 with
/// dft_plus_u=0 must therefore return cleanly.
TEST_F(DFTUBaseTest, AppendSnapshotGuardDftuOffBeforeIterValidation)
{
    Plus_U_Base dftu;
    this->init_dftu(dftu, false);

    // out_occ_mat=true but dft_plus_u=0: the guard must catch dft_plus_u
    // and return before the iter < 1 check fires.
    const DFTU_BASE::OccmatOutputCfg cfg = {1, 1, 5, true, 0};
    const DFTU_BASE::OccmatStepInfo step = {0, 0, 1e-6, 0.5};

    // Use a unique subdirectory so files from other tests in the same
    // binary do not cause false positives.
    const std::string out_dir = testing::TempDir() + "guard_test/";

    // iter=0 is intentionally invalid; the guard-first check ensures
    // this does not abort.
    DFTU_BASE::append_ion_step_snapshot(dftu,
                                        ucell,
                                        out_dir,
                                        step,
                                        false,
                                        false,
                                        cfg,
                                        DFTU_BASE::SOC_LAYOUT_PAULI);

    // If we reach here the guard worked.  No file should have been written
    // (the directory may not even exist, which is fine — the guard
    // returned before any file IO).
    const std::string fn_g1 = out_dir + "occ_matg1.txt";
    std::ifstream ifs(fn_g1.c_str());
    EXPECT_FALSE(ifs.is_open());

    // Sanity check: with dft_plus_u > 0 and iter=0 the function WOULD
    // abort.  We cannot test that path in-process (WARNING_QUIT calls
    // exit/MPI_Abort), so the guard-first test above is the best we can do.
}
