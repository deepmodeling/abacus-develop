#include "gtest/gtest.h"
#include <fstream>
// #include <sstream>
#include <vector>

#define private public
#include "source_io/module_parameter/parameter.h"
#undef private

#ifdef __MPI
#include <mpi.h>
#endif

#define private public
#include "source_cell/sep_cell.h"
#include "source_cell/unitcell.h"
#undef private
pseudo::pseudo()
{
}
pseudo::~pseudo()
{
}
Atom_pseudo::Atom_pseudo()
{
}
Atom_pseudo::~Atom_pseudo()
{
}
Atom::Atom()
{
}
Atom::~Atom()
{
}
InfoNonlocal::InfoNonlocal()
{
}
InfoNonlocal::~InfoNonlocal()
{
}
LCAO_Orbitals::LCAO_Orbitals()
{
}
LCAO_Orbitals::~LCAO_Orbitals()
{
}
Magnetism::Magnetism()
{
}
Magnetism::~Magnetism()
{
}
UnitCell::UnitCell()
{
}
UnitCell::~UnitCell()
{
}

// Test fixture for Sep_Cell tests
class SepCellTest : public ::testing::Test
{
  protected:
    Sep_Cell sep_cell;
    UnitCell ucell;

    // Names for temporary files used in tests
    std::string stru_filename = "STRU_LiF";
    std::string stru_noLi_filename = "STRU_LiF_Warning1";
    std::string f_sep_filename = "F_pbe_50.sep";
    std::string pp_dir = "support/"; // Directory for pseudopotential files

    void SetUp() override
    {
        // Initialize UnitCell for tests that need it.
        // This setup is common for many read_sep_potentials tests.
        ucell.ntype = 2;
        ucell.atom_label.resize(ucell.ntype);
        ucell.atom_label[0] = "Li";
        ucell.atom_label[1] = "F";
        ucell.atoms = new Atom[ucell.ntype];
        ucell.atoms[0].label = "Li";
        ucell.atoms[0].na = 1; // Number of atoms of this type
        ucell.atoms[1].label = "F";
        ucell.atoms[1].na = 1;
    }

    void TearDown() override
    {
        delete[] ucell.atoms;
        ucell.atoms = nullptr;
    }
};

TEST_F(SepCellTest, Constructor)
{
    EXPECT_EQ(sep_cell.get_ntype(), 0);
    EXPECT_DOUBLE_EQ(sep_cell.get_omega(), 0.0);
    EXPECT_DOUBLE_EQ(sep_cell.get_tpiba2(), 0.0);
    EXPECT_TRUE(sep_cell.get_seps().empty());
    EXPECT_TRUE(sep_cell.get_sep_enable().empty());
}

TEST_F(SepCellTest, Init)
{
    sep_cell.init(2);
    EXPECT_EQ(sep_cell.get_ntype(), 2);
    ASSERT_EQ(sep_cell.get_seps().size(), 2);
    ASSERT_EQ(sep_cell.get_sep_enable().size(), 2);
    EXPECT_FALSE(sep_cell.get_sep_enable()[0]);
    EXPECT_FALSE(sep_cell.get_sep_enable()[1]);
    // Check default values of SepPot within seps
    EXPECT_EQ(sep_cell.get_seps()[0].mesh, 0);
    EXPECT_FALSE(sep_cell.get_seps()[0].is_enable);
}

TEST_F(SepCellTest, SetOmega)
{
    sep_cell.set_omega(100.0, 0.25);
    EXPECT_DOUBLE_EQ(sep_cell.get_omega(), 100.0);
    EXPECT_DOUBLE_EQ(sep_cell.get_tpiba2(), 0.25);
}

TEST_F(SepCellTest, ReadSepPotentialsSuccess)
{
#ifdef __MPI
    if (GlobalV::MY_RANK == 0)
    {
#endif

        std::ifstream ifs(pp_dir + stru_filename);
        ASSERT_TRUE(ifs.is_open());

        sep_cell.init(ucell.ntype);
        std::ofstream ofs_running_dummy("dummy_ofs_running.tmp");
        int result = sep_cell.read_sep_potentials(ifs, pp_dir, ofs_running_dummy, ucell.atom_label);
        ifs.close();
        std::remove("dummy_ofs_running.tmp");

        EXPECT_EQ(result, 1); // Expect success (true)

        // Due to the bug mentioned (this->sep_enable[i] is always false),
        // SEP data won't actually be loaded.
        ASSERT_EQ(sep_cell.get_sep_enable().size(), 2);
        EXPECT_FALSE(sep_cell.get_sep_enable()[0]); // Stays false from init
        EXPECT_TRUE(sep_cell.get_sep_enable()[1]);  // Stays false from init

        const auto& seps = sep_cell.get_seps();
        ASSERT_EQ(seps.size(), 2);
        EXPECT_FALSE(seps[0].is_enable); // Default value, not set from file
        EXPECT_EQ(seps[0].mesh, 0);      // Default value
        EXPECT_EQ(seps[0].label, "");    // Default value

        EXPECT_TRUE(seps[1].is_enable); // Default value
        EXPECT_EQ(seps[1].mesh, 1038);  // Default value
        EXPECT_EQ(seps[1].label, "F");  // Default value
        EXPECT_DOUBLE_EQ(seps[1].r_in, 0.0);
        EXPECT_DOUBLE_EQ(seps[1].r_out, 2.5);
        EXPECT_DOUBLE_EQ(seps[1].r_power, 20.0);
        EXPECT_DOUBLE_EQ(seps[1].enhence_a, 1.0);
#ifdef __MPI
    }
    // If run in MPI, other ranks might need to know the outcome or have sep_cell state consistent.
    // For this specific test, only rank 0 performs the read.
    // A broadcast test would cover data consistency across ranks.
#endif
}

TEST_F(SepCellTest, ReadSepPotentialsNoSepFilesSection)
{
#ifdef __MPI
    if (GlobalV::MY_RANK == 0)
    {
#endif

        std::ifstream ifs(pp_dir + stru_noLi_filename);
        ASSERT_TRUE(ifs.is_open());
        std::ofstream ofs_running_dummy("dummy_ofs_running.tmp");

        sep_cell.init(ucell.ntype);
        int result = sep_cell.read_sep_potentials(ifs, pp_dir, ofs_running_dummy, ucell.atom_label);
        ifs.close();
        std::remove("dummy_ofs_running.tmp");

        EXPECT_EQ(result, 0); // Expect failure (false) because "SEP_FILES" not found
#ifdef __MPI
    }
#endif
}

#ifdef __MPI
TEST_F(SepCellTest, BcastSepCell)
{
    sep_cell.init(2); // ntype = 2
    // Rank 0 prepares some data (or reads from file)
    if (GlobalV::MY_RANK == 0)
    {
        sep_cell.set_omega(150.0, 0.75);
        std::ifstream ifs(pp_dir + stru_filename);
        ASSERT_TRUE(ifs.is_open());

        sep_cell.init(ucell.ntype);
        std::ofstream ofs_running_dummy("dummy_ofs_running.tmp");
        int result = sep_cell.read_sep_potentials(ifs, pp_dir, ofs_running_dummy, ucell.atom_label);
        ifs.close();
        std::remove("dummy_ofs_running.tmp");

        EXPECT_EQ(result, 1); // Expect success (true)
    }

    sep_cell.bcast_sep_cell();

    // All ranks should have the same data
    EXPECT_EQ(sep_cell.get_ntype(), 2);
    // Omega and tpiba2 are NOT part of Sep_Cell::bcast_sep_cell, so they remain default on non-zero ranks
    if (GlobalV::MY_RANK == 0)
    {
        EXPECT_DOUBLE_EQ(sep_cell.get_omega(), 150.0);
        EXPECT_DOUBLE_EQ(sep_cell.get_tpiba2(), 0.75);
    }
    else
    {
        EXPECT_DOUBLE_EQ(sep_cell.get_omega(), 0.0);  // Default
        EXPECT_DOUBLE_EQ(sep_cell.get_tpiba2(), 0.0); // Default
    }

    ASSERT_EQ(sep_cell.get_sep_enable().size(), 2);
    // sep_enable will be broadcast as false from rank 0 due to read_sep_potentials bug
    EXPECT_FALSE(sep_cell.get_sep_enable()[0]);
    EXPECT_TRUE(sep_cell.get_sep_enable()[1]);

    const auto& seps = sep_cell.get_seps();
    ASSERT_EQ(seps.size(), 2);

    // Check SepPot data (will be default values due to bug and current test setup)
    EXPECT_EQ(seps[0].label, "");    // Default broadcasted
    EXPECT_EQ(seps[0].mesh, 0);      // Default broadcasted
    EXPECT_FALSE(seps[0].is_enable); // Default broadcasted

    EXPECT_EQ(seps[1].label, "F");  // Default broadcasted
    EXPECT_EQ(seps[1].mesh, 1038);  // Default broadcasted
    EXPECT_TRUE(seps[1].is_enable); // Default broadcasted
    // Note: SepPot::bcast_sep() allocates memory for r and rv on all ranks
    // whenever mesh > 0, regardless of is_enable status
    if (seps[0].mesh > 0)
    {
        EXPECT_NE(seps[0].r, nullptr);
        EXPECT_NE(seps[0].rv, nullptr);
    }
    else
    {
        EXPECT_EQ(seps[0].r, nullptr);
        EXPECT_EQ(seps[0].rv, nullptr);
    }
    EXPECT_NE(seps[1].r, nullptr);
    EXPECT_NE(seps[1].rv, nullptr);
    EXPECT_DOUBLE_EQ(seps[1].r[0], 3.4643182373e-06);
    EXPECT_DOUBLE_EQ(seps[1].rv[0], -2.0868200000e-05);
    EXPECT_DOUBLE_EQ(seps[1].r[7], 2.8965849122e-05);
    EXPECT_DOUBLE_EQ(seps[1].rv[7], -1.9723800000e-05);
}
#endif // __MPI

// Main function for running tests
int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &GlobalV::NPROC);
    MPI_Comm_rank(MPI_COMM_WORLD, &GlobalV::MY_RANK);
#endif

    testing::InitGoogleTest(&argc, argv);

    // Potentially initialize GlobalV::ofs_running here if not handled by test infra
    // e.g., if (GlobalV::MY_RANK == 0) GlobalV::ofs_running.open("sep_cell_test.log");
    // For now, assume it's usable or output to console/dev_null is acceptable.

    int result = RUN_ALL_TESTS();

#ifdef __MPI
    MPI_Finalize();
#endif
    return result;
}
