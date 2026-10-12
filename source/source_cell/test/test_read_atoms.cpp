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
 *  unit test of read_atoms.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - ReadAtomPositionsS1
 *     - read_atom_positions(): spin 1 case
 *   - ReadAtomPositionsS2
 *     - read_atom_positions(): spin 2 case
 *   - ReadAtomPositionsS4Noncolin
 *     - read_atom_positions(): spin 4 noncolinear case
 *   - ReadAtomPositionsS4Colin
 *     - read_atom_positions(): spin 4 colinear case
 *   - ReadAtomPositionsC
 *     - read_atom_positions(): Cartesian coordinates
 *   - ReadAtomPositionsCA
 *     - read_atom_positions(): Cartesian_angstrom coordinates
 *   - ReadAtomPositionsCACXY
 *     - read_atom_positions(): Cartesian_angstrom_center_xy coordinates
 *   - ReadAtomPositionsCACXZ
 *     - read_atom_positions(): Cartesian_angstrom_center_xz coordinates
 *   - ReadAtomPositionsCACYZ
 *     - read_atom_positions(): Cartesian_angstrom_center_yz coordinates
 *   - ReadAtomPositionsCACXYZ
 *     - read_atom_positions(): Cartesian_angstrom_center_xyz coordinates
 *   - ReadAtomPositionsCAU
 *     - read_atom_positions(): Cartesian_au coordinates
 *   - ReadAtomPositionsAutosetMag
 *     - read_atom_positions(): zero-moment start with warning for nspin=4
 *   - ReadAtomPositionsWarning1
 *     - read_atom_positions(): unknown type of coordinates
 *   - ReadAtomPositionsWarning2
 *     - read_atom_positions(): atomic label inconsistency between ATOM_POSITIONS
 *                              and ATOM_SPECIES
 *   - ReadAtomPositionsWarning3
 *     - read_atom_positions(): warning : atom number < 0
 *   - ReadAtomPositionsWarning4
 *     - read_atom_positions(): mismatch in atom number for atom type
 *   - ReadAtomPositionsWarning5
 *     - read_atom_positions(): no atoms can move in MD simulations!
 */

class ReadAtomsTest : public ::testing::Test
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

TEST_F(ReadAtomsTest, ReadAtomPositionsS1)
{
    std::string fn = "./support/STRU_MgO";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    std::ofstream ofs_warning;
    ofs_running.open("read_atom_positions.tmp");
    ofs_warning.open("read_atom_positions.warn");
    // mandatory preliminaries
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    const std::string basis_type = "lcao";
    const std::string orbital_dir = "";
    const std::string init_wfc = "";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = true;
    const bool rpa = false;
    const int nspin = 1;
    const bool fixed_atoms = false;
    const bool noncolin = false;
    const std::string calculation = "scf";
    const std::string esolver_type = "ksdft";
    EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, deepks_setorb, rpa));
    EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat));
    EXPECT_DOUBLE_EQ(ucell->latvec.e11, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e22, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e33, 4.27957);
    unitcell::read_atom_positions(*ucell, ifa, ofs_running, ofs_warning, nspin,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, fixed_atoms, noncolin,
        calculation, esolver_type, 0);
    ofs_running.close();
    ofs_warning.close();
    ifa.close();
    remove("read_atom_positions.tmp");
    remove("read_atom_positions.warn");
}

TEST_F(ReadAtomsTest, ReadAtomPositionsS2)
{
    std::string fn = "./support/STRU_MgO";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    std::ofstream ofs_warning;
    ofs_running.open("read_atom_positions.tmp");
    ofs_warning.open("read_atom_positions.warn");
    // mandatory preliminaries
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    const std::string basis_type = "lcao";
    const std::string orbital_dir = "";
    const std::string init_wfc = "";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = true;
    const bool rpa = false;
    const int nspin = 2;
    const bool fixed_atoms = false;
    const bool noncolin = false;
    const std::string calculation = "scf";
    const std::string esolver_type = "ksdft";
    EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, deepks_setorb, rpa));
    EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat));
    EXPECT_DOUBLE_EQ(ucell->latvec.e11, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e22, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e33, 4.27957);
    unitcell::read_atom_positions(*ucell, ifa, ofs_running, ofs_warning, nspin,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, fixed_atoms, noncolin,
        calculation, esolver_type, 0);
    ofs_running.close();
    ofs_warning.close();
    ifa.close();
    remove("read_atom_positions.tmp");
    remove("read_atom_positions.warn");
}

TEST_F(ReadAtomsTest, ReadAtomPositionsS4Noncolin)
{
    std::string fn = "./support/STRU_MgO";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    std::ofstream ofs_warning;
    ofs_running.open("read_atom_positions.tmp");
    ofs_warning.open("read_atom_positions.warn");
    // mandatory preliminaries
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    const std::string basis_type = "lcao";
    const std::string orbital_dir = "";
    const std::string init_wfc = "";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = true;
    const bool rpa = false;
    const int nspin = 4;
    const bool fixed_atoms = false;
    const bool noncolin = true;
    const std::string calculation = "scf";
    const std::string esolver_type = "ksdft";
    EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, deepks_setorb, rpa));
    EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat));
    EXPECT_DOUBLE_EQ(ucell->latvec.e11, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e22, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e33, 4.27957);
    unitcell::read_atom_positions(*ucell, ifa, ofs_running, ofs_warning, nspin,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, fixed_atoms, noncolin,
        calculation, esolver_type, 0);
    ofs_running.close();
    ofs_warning.close();
    ifa.close();
    remove("read_atom_positions.tmp");
    remove("read_atom_positions.warn");
}

TEST_F(ReadAtomsTest, ReadAtomPositionsS4Colin)
{
    std::string fn = "./support/STRU_MgO";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    std::ofstream ofs_warning;
    ofs_running.open("read_atom_positions.tmp");
    ofs_warning.open("read_atom_positions.warn");
    // mandatory preliminaries
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    const std::string basis_type = "lcao";
    const std::string orbital_dir = "";
    const std::string init_wfc = "";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = true;
    const bool rpa = false;
    const int nspin = 4;
    const bool fixed_atoms = false;
    const bool noncolin = false;
    const std::string calculation = "scf";
    const std::string esolver_type = "ksdft";
    EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, deepks_setorb, rpa));
    EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat));
    EXPECT_DOUBLE_EQ(ucell->latvec.e11, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e22, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e33, 4.27957);
    unitcell::read_atom_positions(*ucell, ifa, ofs_running, ofs_warning, nspin,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, fixed_atoms, noncolin,
        calculation, esolver_type, 0);
    ofs_running.close();
    ofs_warning.close();
    ifa.close();
    remove("read_atom_positions.tmp");
    remove("read_atom_positions.warn");
}

TEST_F(ReadAtomsTest, ReadAtomPositionsC)
{
    std::string fn = "./support/STRU_MgO_c";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    std::ofstream ofs_warning;
    ofs_running.open("read_atom_positions.tmp");
    ofs_warning.open("read_atom_positions.warn");
    // mandatory preliminaries
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    const std::string basis_type = "lcao";
    const std::string orbital_dir = "";
    const std::string init_wfc = "";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = true;
    const bool rpa = false;
    const int nspin = 1;
    const bool fixed_atoms = false;
    const bool noncolin = false;
    const std::string calculation = "scf";
    const std::string esolver_type = "ksdft";
    EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, deepks_setorb, rpa));
    EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat));
    EXPECT_DOUBLE_EQ(ucell->latvec.e11, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e22, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e33, 4.27957);
    unitcell::read_atom_positions(*ucell, ifa, ofs_running, ofs_warning, nspin,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, fixed_atoms, noncolin,
        calculation, esolver_type, 0);
    ofs_running.close();
    ofs_warning.close();
    ifa.close();
    remove("read_atom_positions.tmp");
    remove("read_atom_positions.warn");
}

TEST_F(ReadAtomsTest, ReadAtomPositionsCA)
{
    std::string fn = "./support/STRU_MgO_ca";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    std::ofstream ofs_warning;
    ofs_running.open("read_atom_positions.tmp");
    ofs_warning.open("read_atom_positions.warn");
    // mandatory preliminaries
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    const std::string basis_type = "lcao";
    const std::string orbital_dir = "";
    const std::string init_wfc = "";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = true;
    const bool rpa = false;
    const int nspin = 1;
    const bool fixed_atoms = false;
    const bool noncolin = false;
    const std::string calculation = "scf";
    const std::string esolver_type = "ksdft";
    EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, deepks_setorb, rpa));
    EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat));
    EXPECT_DOUBLE_EQ(ucell->latvec.e11, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e22, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e33, 4.27957);
    unitcell::read_atom_positions(*ucell, ifa, ofs_running, ofs_warning, nspin,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, fixed_atoms, noncolin,
        calculation, esolver_type, 0);
    ofs_running.close();
    ofs_warning.close();
    ifa.close();
    remove("read_atom_positions.tmp");
    remove("read_atom_positions.warn");
}

TEST_F(ReadAtomsTest, ReadAtomPositionsCACXY)
{
    std::string fn = "./support/STRU_MgO_cacxy";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    std::ofstream ofs_warning;
    ofs_running.open("read_atom_positions.tmp");
    ofs_warning.open("read_atom_positions.warn");
    // mandatory preliminaries
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    const std::string basis_type = "lcao";
    const std::string orbital_dir = "";
    const std::string init_wfc = "";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = true;
    const bool rpa = false;
    const int nspin = 1;
    const bool fixed_atoms = false;
    const bool noncolin = false;
    const std::string calculation = "scf";
    const std::string esolver_type = "ksdft";
    EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, deepks_setorb, rpa));
    EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat));
    EXPECT_DOUBLE_EQ(ucell->latvec.e11, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e22, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e33, 4.27957);
    unitcell::read_atom_positions(*ucell, ifa, ofs_running, ofs_warning, nspin,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, fixed_atoms, noncolin,
        calculation, esolver_type, 0);
    ofs_running.close();
    ofs_warning.close();
    ifa.close();
    remove("read_atom_positions.tmp");
    remove("read_atom_positions.warn");
}

TEST_F(ReadAtomsTest, ReadAtomPositionsCACXZ)
{
    std::string fn = "./support/STRU_MgO_cacxz";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    std::ofstream ofs_warning;
    ofs_running.open("read_atom_positions.tmp");
    ofs_warning.open("read_atom_positions.warn");
    // mandatory preliminaries
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    const std::string basis_type = "lcao";
    const std::string orbital_dir = "";
    const std::string init_wfc = "";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = true;
    const bool rpa = false;
    const int nspin = 1;
    const bool fixed_atoms = false;
    const bool noncolin = false;
    const std::string calculation = "scf";
    const std::string esolver_type = "ksdft";
    EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, deepks_setorb, rpa));
    EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat));
    EXPECT_DOUBLE_EQ(ucell->latvec.e11, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e22, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e33, 4.27957);
    unitcell::read_atom_positions(*ucell, ifa, ofs_running, ofs_warning, nspin,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, fixed_atoms, noncolin,
        calculation, esolver_type, 0);
    ofs_running.close();
    ofs_warning.close();
    ifa.close();
    remove("read_atom_positions.tmp");
    remove("read_atom_positions.warn");
}

TEST_F(ReadAtomsTest, ReadAtomPositionsCACYZ)
{
    std::string fn = "./support/STRU_MgO_cacyz";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    std::ofstream ofs_warning;
    ofs_running.open("read_atom_positions.tmp");
    ofs_warning.open("read_atom_positions.warn");
    // mandatory preliminaries
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    const std::string basis_type = "lcao";
    const std::string orbital_dir = "";
    const std::string init_wfc = "";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = true;
    const bool rpa = false;
    const int nspin = 1;
    const bool fixed_atoms = false;
    const bool noncolin = false;
    const std::string calculation = "scf";
    const std::string esolver_type = "ksdft";
    EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, deepks_setorb, rpa));
    EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat));
    EXPECT_DOUBLE_EQ(ucell->latvec.e11, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e22, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e33, 4.27957);
    unitcell::read_atom_positions(*ucell, ifa, ofs_running, ofs_warning, nspin,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, fixed_atoms, noncolin,
        calculation, esolver_type, 0);
    ofs_running.close();
    ofs_warning.close();
    ifa.close();
    remove("read_atom_positions.tmp");
    remove("read_atom_positions.warn");
}

TEST_F(ReadAtomsTest, ReadAtomPositionsCACXYZ)
{
    std::string fn = "./support/STRU_MgO_cacxyz";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    std::ofstream ofs_warning;
    ofs_running.open("read_atom_positions.tmp");
    ofs_warning.open("read_atom_positions.warn");
    // mandatory preliminaries
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    const std::string basis_type = "lcao";
    const std::string orbital_dir = "";
    const std::string init_wfc = "";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = true;
    const bool rpa = false;
    const int nspin = 1;
    const bool fixed_atoms = false;
    const bool noncolin = false;
    const std::string calculation = "scf";
    const std::string esolver_type = "ksdft";
    EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, deepks_setorb, rpa));
    EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat));
    EXPECT_DOUBLE_EQ(ucell->latvec.e11, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e22, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e33, 4.27957);
    unitcell::read_atom_positions(*ucell, ifa, ofs_running, ofs_warning, nspin,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, fixed_atoms, noncolin,
        calculation, esolver_type, 0);
    ofs_running.close();
    ofs_warning.close();
    ifa.close();
    remove("read_atom_positions.tmp");
    remove("read_atom_positions.warn");
}

TEST_F(ReadAtomsTest, ReadAtomPositionsCAU)
{
    std::string fn = "./support/STRU_MgO_cau";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    std::ofstream ofs_warning;
    ofs_running.open("read_atom_positions.tmp");
    ofs_warning.open("read_atom_positions.warn");
    // mandatory preliminaries
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    const std::string basis_type = "lcao";
    const std::string orbital_dir = "";
    const std::string init_wfc = "";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = true;
    const bool rpa = false;
    const int nspin = 1;
    const bool fixed_atoms = true;
    const bool noncolin = false;
    const std::string calculation = "scf";
    const std::string esolver_type = "ksdft";
    EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, deepks_setorb, rpa));
    EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat));
    EXPECT_DOUBLE_EQ(ucell->latvec.e11, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e22, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e33, 4.27957);
    unitcell::read_atom_positions(*ucell, ifa, ofs_running, ofs_warning, nspin,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, fixed_atoms, noncolin,
        calculation, esolver_type, 0);
    ofs_running.close();
    ofs_warning.close();
    ifa.close();
    remove("read_atom_positions.tmp");
    remove("read_atom_positions.warn");
}

TEST_F(ReadAtomsTest, ReadAtomPositionsAutosetMag)
{
    std::string fn = "./support/STRU_MgO";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    std::ofstream ofs_warning;
    ofs_running.open("read_atom_positions.tmp");
    ofs_warning.open("read_atom_positions.warn");
    // mandatory preliminaries
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    const std::string basis_type = "lcao";
    const std::string orbital_dir = "";
    const std::string init_wfc = "";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = true;
    const bool rpa = false;
    const bool fixed_atoms = false;
    const bool noncolin = false;
    const std::string calculation = "scf";
    const std::string esolver_type = "ksdft";
    int nspin = 2;
    EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, deepks_setorb, rpa));
    EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat));
    EXPECT_DOUBLE_EQ(ucell->latvec.e11, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e22, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e33, 4.27957);
    unitcell::read_atom_positions(*ucell, ifa, ofs_running, ofs_warning, nspin,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, fixed_atoms, noncolin,
        calculation, esolver_type, 0);
    for (int it = 0; it < ucell->ntype; it++)
    {
        for (int ia = 0; ia < ucell->atoms[it].na; ia++)
        {
            EXPECT_DOUBLE_EQ(ucell->atoms[it].mag[ia], 1.0);
            EXPECT_DOUBLE_EQ(ucell->atoms[it].m_loc_[ia].z, 1.0);
        }
    }
    // for nspin == 4
    // Issue #5939: nspin=4 with no mag in STRU no longer autosets (1,1,1);
    // all moments stay zero and a warning is emitted instead.
    nspin = 4;
    testing::internal::CaptureStdout();
    unitcell::read_atom_positions(*ucell, ifa, ofs_running, ofs_warning, nspin,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, fixed_atoms, noncolin,
        calculation, esolver_type, 0);
    const std::string stdout_output = testing::internal::GetCapturedStdout();
    for (int it = 0; it < ucell->ntype; it++)
    {
        for (int ia = 0; ia < ucell->atoms[it].na; ia++)
        {
            EXPECT_DOUBLE_EQ(ucell->atoms[it].mag[ia], 0.0);
            EXPECT_DOUBLE_EQ(ucell->atoms[it].m_loc_[ia].x, 0.0);
            EXPECT_DOUBLE_EQ(ucell->atoms[it].m_loc_[ia].y, 0.0);
            EXPECT_DOUBLE_EQ(ucell->atoms[it].m_loc_[ia].z, 0.0);
        }
    }
    // The zero-moment warning must reach both stdout and the running log.
    EXPECT_NE(stdout_output.find("no initial magnetization is set in STRU"),
              std::string::npos);
    ofs_running.flush();
    std::ifstream ifs_log("read_atom_positions.tmp");
    std::string log_content;
    std::string log_line;
    while (std::getline(ifs_log, log_line))
    {
        log_content += log_line;
    }
    ifs_log.close();
    EXPECT_NE(log_content.find("no initial magnetization is set in STRU"),
              std::string::npos);
    ofs_running.close();
    ofs_warning.close();
    ifa.close();
    remove("read_atom_positions.tmp");
    remove("read_atom_positions.warn");
}

TEST_F(ReadAtomsTest, ReadAtomPositionsWarning1)
{
    std::string fn = "./support/STRU_MgO_WarningC1";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    std::ofstream ofs_warning;
    ofs_running.open("read_atom_positions.tmp");
    ofs_warning.open("read_atom_positions.warn");
    // mandatory preliminaries
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    const std::string basis_type = "lcao";
    const std::string orbital_dir = "";
    const std::string init_wfc = "";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = true;
    const bool rpa = false;
    const int nspin = 1;
    const bool fixed_atoms = false;
    const bool noncolin = false;
    const std::string calculation = "scf";
    const std::string esolver_type = "ksdft";
    EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, deepks_setorb, rpa));
    EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat));
    EXPECT_DOUBLE_EQ(ucell->latvec.e11, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e22, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e33, 4.27957);
    EXPECT_NO_THROW(unitcell::read_atom_positions(*ucell, ifa, ofs_running, ofs_warning, nspin,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, fixed_atoms, noncolin,
        calculation, esolver_type, 0));
    ofs_running.close();
    ofs_warning.close();
    ifa.close();
    // check warning file
    std::ifstream ifs_tmp;
    ifs_tmp.open("read_atom_positions.warn");
    std::string str((std::istreambuf_iterator<char>(ifs_tmp)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("There are several options for you:"));
    EXPECT_THAT(str, testing::HasSubstr("Direct"));
    EXPECT_THAT(str, testing::HasSubstr("Cartesian_angstrom"));
    EXPECT_THAT(str, testing::HasSubstr("Cartesian_au"));
    EXPECT_THAT(str, testing::HasSubstr("Cartesian_angstrom_center_xy"));
    EXPECT_THAT(str, testing::HasSubstr("Cartesian_angstrom_center_xz"));
    EXPECT_THAT(str, testing::HasSubstr("Cartesian_angstrom_center_yz"));
    EXPECT_THAT(str, testing::HasSubstr("Cartesian_angstrom_center_xyz"));
    ifs_tmp.close();
    remove("read_atom_positions.tmp");
    remove("read_atom_positions.warn");
}

TEST_F(ReadAtomsTest, ReadAtomPositionsWarning2)
{
    std::string fn = "./support/STRU_MgO_WarningC2";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    std::ofstream ofs_warning;
    ofs_running.open("read_atom_positions.tmp");
    ofs_warning.open("read_atom_positions.warn");
    // mandatory preliminaries
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    const std::string basis_type = "lcao";
    const std::string orbital_dir = "";
    const std::string init_wfc = "";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = true;
    const bool rpa = false;
    const int nspin = 1;
    const bool fixed_atoms = false;
    const bool noncolin = false;
    const std::string calculation = "scf";
    const std::string esolver_type = "ksdft";
    EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, deepks_setorb, rpa));
    EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat));
    EXPECT_DOUBLE_EQ(ucell->latvec.e11, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e22, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e33, 4.27957);
    EXPECT_NO_THROW(unitcell::read_atom_positions(*ucell, ifa, ofs_running, ofs_warning, nspin,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, fixed_atoms, noncolin,
        calculation, esolver_type, 0));
    ofs_running.close();
    ofs_warning.close();
    ifa.close();
    // check warning file
    std::ifstream ifs_tmp;
    ifs_tmp.open("read_atom_positions.warn");
    std::string str((std::istreambuf_iterator<char>(ifs_tmp)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("Label read from ATOMIC_POSITIONS is Mo"));
    EXPECT_THAT(str, testing::HasSubstr("Label from ATOMIC_SPECIES is Mg"));
    ifs_tmp.close();
    remove("read_atom_positions.tmp");
    remove("read_atom_positions.warn");
}

TEST_F(ReadAtomsTest, ReadAtomPositionsWarning3)
{
    std::string fn = "./support/STRU_MgO_WarningC3";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    ofs_running.open("read_atom_positions.tmp");
    GlobalV::ofs_warning.open("read_atom_positions.warn");
    // mandatory preliminaries
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    const std::string basis_type = "lcao";
    const std::string orbital_dir = "";
    const std::string init_wfc = "";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = true;
    const bool rpa = false;
    const int nspin = 1;
    const bool fixed_atoms = false;
    const bool noncolin = false;
    const std::string calculation = "scf";
    const std::string esolver_type = "ksdft";
    EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, deepks_setorb, rpa));
    EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat));
    EXPECT_DOUBLE_EQ(ucell->latvec.e11, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e22, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e33, 4.27957);
    EXPECT_NO_THROW(unitcell::read_atom_positions(*ucell, ifa, ofs_running, GlobalV::ofs_warning, nspin,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, fixed_atoms, noncolin,
        calculation, esolver_type, 0));
    ofs_running.close();
    GlobalV::ofs_warning.close();
    ifa.close();
    // check warning file
    std::ifstream ifs_tmp;
    ifs_tmp.open("read_atom_positions.warn");
    std::string str((std::istreambuf_iterator<char>(ifs_tmp)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("read_atom_positions  warning :  atom number < 0."));
    ifs_tmp.close();
    remove("read_atom_positions.tmp");
    remove("read_atom_positions.warn");
}

TEST_F(ReadAtomsTest, ReadAtomPositionsWarning4)
{
    std::string fn = "./support/STRU_MgO_WarningC4";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    std::ofstream ofs_warning;
    ofs_running.open("read_atom_positions.tmp");
    ofs_warning.open("read_atom_positions.warn");
    // mandatory preliminaries
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->orbital_fn.resize(ucell->ntype);
    ucell->set_atom_flag = true;
    const std::string basis_type = "lcao";
    const std::string orbital_dir = "";
    const std::string init_wfc = "";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = true;
    const bool rpa = false;
    const int nspin = 1;
    const bool fixed_atoms = false;
    const bool noncolin = false;
    const std::string calculation = "scf";
    const std::string esolver_type = "ksdft";
    EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, deepks_setorb, rpa));
    EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat));
    EXPECT_DOUBLE_EQ(ucell->latvec.e11, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e22, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e33, 4.27957);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(unitcell::read_atom_positions(*ucell, ifa, ofs_running, ofs_warning, nspin,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, fixed_atoms, noncolin,
        calculation, esolver_type, 0), ::testing::ExitedWithCode(1), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("read_atom_positions, mismatch in atom number for atom type: Mg"));
    ofs_running.close();
    ofs_warning.close();
    ifa.close();
    remove("read_atom_positions.tmp");
    remove("read_atom_positions.warn");
}

TEST_F(ReadAtomsTest, ReadAtomPositionsWarning5)
{
    std::string fn = "./support/STRU_MgO";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    ofs_running.open("read_atom_positions.tmp");
    GlobalV::ofs_warning.open("read_atom_positions.warn");
    // mandatory preliminaries
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    const std::string basis_type = "lcao";
    const std::string orbital_dir = "";
    const std::string init_wfc = "";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = true;
    const bool rpa = false;
    const int nspin = 1;
    const bool fixed_atoms = true;
    const bool noncolin = false;
    const std::string calculation = "md";
    const std::string esolver_type = "ksdft";
    EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, deepks_setorb, rpa));
    EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat));
    EXPECT_DOUBLE_EQ(ucell->latvec.e11, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e22, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e33, 4.27957);
    EXPECT_NO_THROW(unitcell::read_atom_positions(*ucell, ifa, ofs_running, GlobalV::ofs_warning, nspin,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, fixed_atoms, noncolin,
        calculation, esolver_type, 0));
    ofs_running.close();
    GlobalV::ofs_warning.close();
    ifa.close();
    // check warning file
    std::ifstream ifs_tmp;
    ifs_tmp.open("read_atom_positions.warn");
    std::string str((std::istreambuf_iterator<char>(ifs_tmp)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("read_atoms  warning : no atoms can move in MD simulations!"));
    ifs_tmp.close();
    remove("read_atom_positions.tmp");
    remove("read_atom_positions.warn");
}

#endif // __LCAO
