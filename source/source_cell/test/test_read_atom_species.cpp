#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "source_base/global_variable.h"
#include "source_base/mathzone.h"
#include "source_cell/unitcell.h"
#include "source_cell/read_stru.h"

// The test-only cell_info object library does not contain magnetism.cpp,
// so the Magnetism constructor/destructor must be provided locally.
Magnetism::Magnetism()
{
    this->tot_mag = 0.0;
    this->abs_mag = 0.0;
}
Magnetism::~Magnetism()
{
}

#ifdef __LCAO
/************************************************
 *  unit test of read_atom_species.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - ReadAtomSpecies
 *     - read_atom_species(): a successful case
 *   - ReadAtomSpeciesWarning1
 *     - read_atom_species(): unrecognized pseudopotential type
 *   - ReadAtomSpeciesLatName
 *     - read_lattice_constant(): all supported lattice names
 *   - ReadLatticeConstantWarning1
 *     - read_lattice_constant(): lattice constant <= 0
 *   - ReadLatticeConstantWarning2
 *     - read_lattice_constant(): LATTICE_PARAMETERS without lattice type
 *   - ReadLatticeConstantWarning3
 *     - read_lattice_constant(): LATTICE_VECTORS with explicit lattice type
 *   - ReadAtomSpeciesWarning5
 *     - read_lattice_constant(): unsupported lattice name
 */

class ReadAtomSpeciesTest : public ::testing::Test
{
  protected:
    std::unique_ptr<UnitCell> ucell{new UnitCell};
    std::string output;
    void SetUp() override
    {
        ucell->ntype = 2;
        ucell->pseudo_fn.resize(ucell->ntype);
        ucell->pseudo_type.resize(ucell->ntype);
        ucell->orbital_fn.resize(ucell->ntype);
    }
    void TearDown() override
    {
        ucell->orbital_fn.shrink_to_fit();
    }
};

using ReadAtomSpeciesDeathTest = ReadAtomSpeciesTest;

TEST_F(ReadAtomSpeciesTest, ReadAtomSpecies)
{
    std::string fn = "./support/STRU_MgO";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    ofs_running.open("read_atom_species.tmp");
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    const std::string basis_type = "lcao";
    const std::string orbital_dir = "";
    const std::string init_wfc = "";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = true;
    const bool rpa = false;
    EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, deepks_setorb, rpa));
    EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running, ucell->lat));
    EXPECT_DOUBLE_EQ(ucell->latvec.e11, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e22, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e33, 4.27957);
    ofs_running.close();
    ifa.close();
    remove("read_atom_species.tmp");
}

TEST_F(ReadAtomSpeciesTest, ReadAtomSpeciesWarning1)
{
    std::string fn = "./support/STRU_MgO_Warning1";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    ofs_running.open("read_atom_species.txt");
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    const std::string basis_type = "lcao";
    const std::string orbital_dir = "";
    const std::string init_wfc = "";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = true;
    const bool rpa = false;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, deepks_setorb, rpa), ::testing::ExitedWithCode(1), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("unrecognized pseudopotential type."));
    ofs_running.close();
    ifa.close();
    //remove("read_atom_species.txt");
}

TEST_F(ReadAtomSpeciesTest, ReadLatticeConstantWarning1)
{
    std::string fn = "./support/STRU_MgO_Warning2";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    ofs_running.open("read_atom_species1.tmp");
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat), ::testing::ExitedWithCode(1), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("Lattice constant <= 0.0"));
    ofs_running.close();
    ifa.close();
    remove("read_atom_species1.tmp");
}

TEST_F(ReadAtomSpeciesTest, ReadLatticeConstantWarning2)
{
    std::string fn = "./support/STRU_MgO_Warning3";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    ofs_running.open("read_atom_species.tmp");
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat), ::testing::ExitedWithCode(1), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output,
                testing::HasSubstr("do not use LATTICE_PARAMETERS without explicit specification of lattice type"));
    ofs_running.close();
    ifa.close();
    remove("read_atom_species.tmp");
}

TEST_F(ReadAtomSpeciesTest, ReadLatticeConstantWarning3)
{
    std::string fn = "./support/STRU_MgO_Warning4";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    ofs_running.open("read_atom_species.tmp");
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    ucell->latName = "bcc";
    testing::internal::CaptureStdout();
    EXPECT_EXIT(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat), ::testing::ExitedWithCode(1), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output,
                testing::HasSubstr("do not use LATTICE_VECTORS along with explicit specification of lattice type"));
    ofs_running.close();
    ifa.close();
    remove("read_atom_species.tmp");
}

TEST_F(ReadAtomSpeciesTest, ReadAtomSpeciesLatName)
{
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    std::vector<std::string> latName_in = {"sc",
                                           "fcc",
                                           "bcc",
                                           "hexagonal",
                                           "trigonal",
                                           "st",
                                           "bct",
                                           "so",
                                           "baco",
                                           "fco",
                                           "bco",
                                           "sm",
                                           "bacm",
                                           "triclinic"};
    for (int i = 0; i < latName_in.size(); ++i)
    {
        std::string fn = "./support/STRU_MgO_LatName";
        std::ifstream ifa(fn.c_str());
        std::ofstream ofs_running;
        ofs_running.open("read_atom_species.tmp");
        ucell->latName = latName_in[i];
        EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat));
        if (ucell->latName == "sc")
        {
            EXPECT_DOUBLE_EQ(ucell->latvec.e11, 1.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e22, 1.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e33, 1.0);
        }
        ofs_running.close();
        ifa.close();
        remove("read_atom_species.tmp");
    }
}

TEST_F(ReadAtomSpeciesDeathTest, ReadAtomSpeciesWarning5)
{
    std::string fn = "./support/STRU_MgO_LatName";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    ofs_running.open("read_atom_species.tmp");
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    ucell->latName = "arbitrary";
    testing::internal::CaptureStdout();
    EXPECT_EXIT(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat), ::testing::ExitedWithCode(1), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("latname not supported"));
    ofs_running.close();
    ifa.close();
    remove("read_atom_species.tmp");
}

#endif // __LCAO
