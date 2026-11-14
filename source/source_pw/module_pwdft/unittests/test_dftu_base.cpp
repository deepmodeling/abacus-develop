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



void DFTUBaseTest::SetUp()
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

void DFTUBaseTest::TearDown()
    {
        // set_atom_flag is false, so ~UnitCell() skips atoms but frees
        // iat2it / iat2ia itself; only atoms must be deleted here.
        delete[] ucell.atoms;
    }

    /// Call init_base with the given Yukawa switch on a fresh d orbital
void DFTUBaseTest::init_dftu(Plus_U_Base& dftu, const bool yukawa_potential)
    {
        init_dftu_for_spin(dftu,
                           yukawa_potential,
                           2,
                           1,
                           DFTU_BASE::SOC_LAYOUT_PAULI);
    }

void DFTUBaseTest::init_dftu_for_spin(Plus_U_Base& dftu,
                            const bool yukawa_potential,
                            const int nspin,
                            const int npol,
                            const DFTU_BASE::OccmatSocLayout& soc_layout)
    {
        const std::vector<int> l_channel = {2};
        const std::vector<double> hubbard_u = {0.0};
        dftu.init_base(ucell,
                       npol,
                       nspin,
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
                       soc_layout);
    }

/// Enabling Yukawa screening constructs its run-lifetime state.
TEST_F(DFTUBaseTest, InitBaseYukawaEnabledCreatesObject)
{
    Plus_U_Base dftu;

    init_dftu(dftu, true);
    EXPECT_TRUE(dftu.use_yukawa());
}

/// Disabling Yukawa screening leaves its optional state empty.
TEST_F(DFTUBaseTest, InitBaseYukawaDisabledLeavesObjectEmpty)
{
    Plus_U_Base dftu;

    init_dftu(dftu, false);
    EXPECT_FALSE(dftu.use_yukawa());
}

/// U-ramping starts from zero and advances toward the configured target.
TEST_F(DFTUBaseTest, InitBaseSeedsURampingState)
{
    Plus_U_Base dftu;
    const std::vector<int> l_channel = {2};
    const std::vector<double> hubbard_u = {5.0};

    dftu.init_base(ucell,
                   1,
                   2,
                   l_channel,
                   false,
                   0.5,
                   "",
                   "",
                   "none",
                   "cpu",
                   hubbard_u,
                   0.5,
                   0,
                   0,
                   DFTU_BASE::SOC_LAYOUT_PAULI);

    EXPECT_DOUBLE_EQ(dftu.get_u_current(0), 0.0);
    dftu.uramping_update();
    EXPECT_DOUBLE_EQ(dftu.get_u_current(0), 0.5);
    EXPECT_FALSE(dftu.u_converged());
    dftu.uramping_update();
    EXPECT_DOUBLE_EQ(dftu.get_u_current(0), 1.0);
    EXPECT_FALSE(dftu.u_converged());
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


/// For init_occ_mat=1 the file-loaded occupation matrix must seed the
/// Hamiltonian at the first SCF step.  Before the fix the whole
/// iter_init_dftu_pw was skipped at the first step, so pot_onsite and
/// E_plusU stayed zero and the file matrix was overwritten by the
/// psi-reprojected matrix at iter=2 without ever entering the Hamiltonian.
/// This test verifies the mechanism: after init_base reads the file,
/// compute_pot_uterm_and_energy must produce a non-zero potential and
/// energy from the file-loaded matrix.
TEST_F(DFTUBaseTest, InitOccMat1BuildsPotentialFromFile)
{
    const std::string dir = testing::TempDir();
    const std::string fn = dir + "occ_mat.txt";

    // One Fe atom, L=2, nspin=2.  All 25 elements per spin set to a
    // constant so expected pot_onsite and energy are easy to compute.
    const double occ_up = 0.1;
    const double occ_dn = 0.2;
    {
        std::ofstream ofs(fn);
        ASSERT_TRUE(ofs.is_open());
        ofs << "# compact test fixture\n";
        ofs << " Fe Atom 1 L 2 mag 0.0\n";
        ofs << " spin 1 nelec 0.5\n";
        for (int m0 = 0; m0 < 5; ++m0)
        {
            for (int m1 = 0; m1 < 5; ++m1) { ofs << " " << occ_up; }
            ofs << "\n";
        }
        ofs << " spin 2 nelec 0.5\n";
        for (int m0 = 0; m0 < 5; ++m0)
        {
            for (int m1 = 0; m1 < 5; ++m1) { ofs << " " << occ_dn; }
            ofs << "\n";
        }
    }

    Plus_U_Base dftu;
    const std::vector<int> l_channel = {2};
    const std::vector<double> hubbard_u = {4.0};
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
                   1,                // init_occ_mat
                   0,                // mixing_dftu
                   DFTU_BASE::SOC_LAYOUT_PAULI);

    ASSERT_TRUE(dftu.is_occmat_ready());
    EXPECT_NEAR(dftu.occmat().get(0, 2, 0, 0, 0), occ_up, 1e-12);
    EXPECT_NEAR(dftu.occmat().get(0, 2, 1, 4, 4), occ_dn, 1e-12);

    // Before building the potential, uterm_mat is all zeros.
    const std::vector<std::complex<double>>& uterm0 = dftu.get_uterm_mat();
    for (size_t i = 0; i < uterm0.size(); ++i)
    {
        EXPECT_EQ(uterm0[i].real(), 0.0);
    }

    // Build pot_onsite and E_plusU from the file-loaded matrix — this is
    // the call the fix adds at the first SCF step for init_occ_mat=1.
    pw::compute_pot_uterm_and_energy(ucell,
                                    2,                  // nspin
                                    dftu.get_u_current_vec(),
                                    dftu.get_l_channel_vec(),
                                    dftu.get_uterm_mat_index(),
                                    dftu.occmat(),
                                    dftu.get_uterm_mat(),
                                    dftu.energy_ref());

    const double u_value = 4.0;
    const int m_size = 5;
    const int size = m_size * m_size;          // 25
    const std::vector<std::complex<double>>& uterm = dftu.get_uterm_mat();
    ASSERT_EQ(uterm.size(), static_cast<size_t>(2 * size));

    // spin-up block [0 .. size): diagonal = U*(0.5 - occ_up) = 1.6
    // off-diagonal = U*(0 - occ_up) = -0.4
    for (int m = 0; m < m_size; ++m)
    {
        EXPECT_NEAR(uterm[m * m_size + m].real(), u_value * (0.5 - occ_up), 1e-12)
            << "spin-up diagonal mismatch at m=" << m;
    }
    EXPECT_NEAR(uterm[0 * m_size + 1].real(), -u_value * occ_up, 1e-12);
    EXPECT_NEAR(uterm[1 * m_size + 0].real(), -u_value * occ_up, 1e-12);

    // spin-down block [size .. 2*size): diagonal = U*(0.5 - occ_dn) = 1.2
    // off-diagonal = U*(0 - occ_dn) = -0.8
    for (int m = 0; m < m_size; ++m)
    {
        EXPECT_NEAR(uterm[size + m * m_size + m].real(),
                    u_value * (0.5 - occ_dn), 1e-12)
            << "spin-down diagonal mismatch at m=" << m;
    }
    EXPECT_NEAR(uterm[size + 0 * m_size + 1].real(), -u_value * occ_dn, 1e-12);

    // E_U = U * 0.5 * 25 * (occ_up^2 + occ_dn^2)
    //     = 4.0 * 0.5 * 25 * (0.01 + 0.04) = 2.5
    const double expected_energy = u_value * 0.5 * size
                                   * (occ_up * occ_up + occ_dn * occ_dn);
    EXPECT_NEAR(dftu.get_energy(), expected_energy, 1e-12);
    EXPECT_GT(dftu.get_energy(), 0.0);
}

/// For init_occ_mat=2 the occupation matrix is fixed from the file
/// throughout the calculation.  Every SCF step must rebuild the potential
/// from the same matrix.  This test verifies that (1) the file matrix
/// produces a non-zero potential and energy, and (2) a second call with
/// the unchanged matrix reproduces the same results — the defining
/// property of mode 2 ("fixed matrix").
TEST_F(DFTUBaseTest, InitOccMat2BuildsPotentialFromFileAndIsStable)
{
    const std::string dir = testing::TempDir();
    const std::string fn = dir + "occ_mat.txt";

    const double occ_up = 0.1;
    const double occ_dn = 0.2;
    {
        std::ofstream ofs(fn);
        ASSERT_TRUE(ofs.is_open());
        ofs << "# compact test fixture\n";
        ofs << " Fe Atom 1 L 2 mag 0.0\n";
        ofs << " spin 1 nelec 0.5\n";
        for (int m0 = 0; m0 < 5; ++m0)
        {
            for (int m1 = 0; m1 < 5; ++m1) { ofs << " " << occ_up; }
            ofs << "\n";
        }
        ofs << " spin 2 nelec 0.5\n";
        for (int m0 = 0; m0 < 5; ++m0)
        {
            for (int m1 = 0; m1 < 5; ++m1) { ofs << " " << occ_dn; }
            ofs << "\n";
        }
    }

    Plus_U_Base dftu;
    const std::vector<int> l_channel = {2};
    const std::vector<double> hubbard_u = {4.0};
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
                   2,                // init_occ_mat — mode 2: fixed matrix
                   0,                // mixing_dftu
                   DFTU_BASE::SOC_LAYOUT_PAULI);

    ASSERT_TRUE(dftu.is_occmat_ready());

    // First call — build potential from the file-loaded matrix.
    pw::compute_pot_uterm_and_energy(ucell,
                                    2,
                                    dftu.get_u_current_vec(),
                                    dftu.get_l_channel_vec(),
                                    dftu.get_uterm_mat_index(),
                                    dftu.occmat(),
                                    dftu.get_uterm_mat(),
                                    dftu.energy_ref());

    const double energy_first = dftu.get_energy();
    EXPECT_GT(energy_first, 0.0);

    // Snapshot the potential after the first build.
    const std::vector<std::complex<double>> uterm_first = dftu.get_uterm_mat();

    // Second call — mode 2 means the matrix never changes, so the
    // potential and energy must be identical.
    pw::compute_pot_uterm_and_energy(ucell,
                                    2,
                                    dftu.get_u_current_vec(),
                                    dftu.get_l_channel_vec(),
                                    dftu.get_uterm_mat_index(),
                                    dftu.occmat(),
                                    dftu.get_uterm_mat(),
                                    dftu.energy_ref());

    const double energy_second = dftu.get_energy();
    EXPECT_NEAR(energy_first, energy_second, 1e-12);

    const std::vector<std::complex<double>>& uterm_second = dftu.get_uterm_mat();
    ASSERT_EQ(uterm_first.size(), uterm_second.size());
    for (size_t i = 0; i < uterm_first.size(); ++i)
    {
        EXPECT_NEAR(uterm_first[i].real(), uterm_second[i].real(), 1e-12)
            << "uterm mismatch at index " << i;
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
