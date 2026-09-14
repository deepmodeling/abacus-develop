#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "memory"
#include "source_base/mathzone.h"
#include "source_base/global_variable.h"
#include "source_cell/unitcell.h"
#include<vector>
#include<valarray>
#include <streambuf>
#include "prepare_unitcell.h"
#include "source_cell/update_cell.h"

Magnetism::Magnetism()
{
    this->tot_mag = 0.0;
    this->abs_mag = 0.0;
}
Magnetism::~Magnetism()
{
}

/************************************************
 *  unit test of class UnitCell
 ***********************************************/

/**
 * - Tested Functions:
 *   - SetupCellS1
 *     - setup_cell: spin 1 case
 *   - SetupCellS2
 *     - setup_cell: spin 2 case
 *   - SetupCellS4
 *     - setup_cell: spin 4 case
 *   - SetupCellWarning1
 *     - setup_cell: Can not find the file containing atom positions.
 *   - SetupCellWarning2
 *     - setup_cell: Something wrong during read_atom_positions
 *   - SetupCellAfterVC
 *     - setup_cell_after_vc
 */

class UcellTest : public ::testing::Test
{
protected:
    std::unique_ptr<UnitCell> ucell{new UnitCell};
    std::string output;

    const double symmetry_prec = 1e-5;
    const int dfthalf_type = 0;
    const std::string pseudo_dir = "./support";
    const std::string basis_type = "pw";
    const std::string orbital_dir = "./";
    const std::string init_wfc = "atomic";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = false;
    const bool rpa = false;
    const bool fixed_atoms = false;
    const bool noncolin = false;
    const std::string calculation = "scf";
    const std::string esolver_type = "cg";

    void SetUp()
    {
        ucell->lmaxmax = 2;
        ucell->ntype   = 2;
        ucell->pseudo_fn.resize(ucell->ntype);
        ucell->pseudo_type.resize(ucell->ntype);
        ucell->orbital_fn.resize(ucell->ntype);
    }
};

using UcellDeathTest = UcellTest;

TEST_F(UcellTest,SetupCellS1)
{
    std::string fn = "./support/STRU_MgO";
    std::ofstream ofs_running;
    ofs_running.open("setup_cell.tmp");
    const int nspin = 1;
    
    ucell->setup_cell(fn, ofs_running, symmetry_prec, dfthalf_type, pseudo_dir, nspin,
        basis_type, orbital_dir, init_wfc, onsite_radius, deepks_setorb, rpa,
        fixed_atoms, noncolin, calculation, esolver_type, 0);
    ofs_running.close();
    remove("setup_cell.tmp");
}

TEST_F(UcellTest,SetupCellS2)
{
    std::string fn = "./support/STRU_MgO";
    std::ofstream ofs_running;
    ofs_running.open("setup_cell.tmp");
    const int nspin = 2;
    
    ucell->setup_cell(fn, ofs_running, symmetry_prec, dfthalf_type, pseudo_dir, nspin,
        basis_type, orbital_dir, init_wfc, onsite_radius, deepks_setorb, rpa,
        fixed_atoms, noncolin, calculation, esolver_type, 0);
    ofs_running.close();
    remove("setup_cell.tmp");
}

TEST_F(UcellTest,SetupCellS4)
{
    std::string fn = "./support/STRU_MgO";
    std::ofstream ofs_running;
    ofs_running.open("setup_cell.tmp");
    const int nspin = 4;
    
    ucell->setup_cell(fn, ofs_running, symmetry_prec, dfthalf_type, pseudo_dir, nspin,
        basis_type, orbital_dir, init_wfc, onsite_radius, deepks_setorb, rpa,
        fixed_atoms, noncolin, calculation, esolver_type, 0);
    ofs_running.close();
    remove("setup_cell.tmp");
}

TEST_F(UcellDeathTest,SetupCellWarning1)
{
    std::string fn = "./STRU_MgO";
    std::ofstream ofs_running;
    ofs_running.open("setup_cell.tmp");
    
    testing::internal::CaptureStdout();
    const int nspin = 1;
    EXPECT_EXIT(ucell->setup_cell(fn, ofs_running, symmetry_prec, dfthalf_type, pseudo_dir, nspin,
        basis_type, orbital_dir, init_wfc, onsite_radius, deepks_setorb, rpa,
        fixed_atoms, noncolin, calculation, esolver_type, 0), ::testing::ExitedWithCode(1), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output,testing::HasSubstr("Can not find the file containing atom positions.!"));
    ofs_running.close();
    remove("setup_cell.tmp");
}

TEST_F(UcellDeathTest,SetupCellWarning2)
{
    std::string fn = "./support/STRU_MgO_WarningC2";
    std::ofstream ofs_running;
    ofs_running.open("setup_cell.tmp");
    
    testing::internal::CaptureStdout();
    const int nspin = 1;
    EXPECT_EXIT(ucell->setup_cell(fn, ofs_running, symmetry_prec, dfthalf_type, pseudo_dir, nspin,
        basis_type, orbital_dir, init_wfc, onsite_radius, deepks_setorb, rpa,
        fixed_atoms, noncolin, calculation, esolver_type, 0), ::testing::ExitedWithCode(1), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output,testing::HasSubstr("Something wrong during read_atom_positions"));
    ofs_running.close();
    remove("setup_cell.tmp");
}

TEST_F(UcellTest,SetupCellAfterVC)
{
    std::string fn = "./support/STRU_MgO";
    std::ofstream ofs_running;
    ofs_running.open("setup_cell.tmp");
    const int nspin = 1;

    ucell->setup_cell(fn, ofs_running, symmetry_prec, dfthalf_type, pseudo_dir, nspin,
        basis_type, orbital_dir, init_wfc, onsite_radius, deepks_setorb, rpa,
        fixed_atoms, noncolin, calculation, esolver_type, 0);
    ucell->lat0 = 1.0;
    ucell->latvec.Zero();
    ucell->latvec.e11 = 10.0;
    ucell->latvec.e22 = 10.0;
    ucell->latvec.e33 = 10.0;
    for (int i =0;i<ucell->ntype;i++)
    {
        ucell->atoms[i].na = 1;
        ucell->atoms[i].taud.resize(ucell->atoms[i].na);
        ucell->atoms[i].tau.resize(ucell->atoms[i].na);
        ucell->atoms[i].taud[0].x = 0.1;
        ucell->atoms[i].taud[0].y = 0.1;
        ucell->atoms[i].taud[0].z = 0.1;
    }
    
    unitcell::setup_cell_after_vc(*ucell,ofs_running, nspin);
    EXPECT_EQ(ucell->lat0_angstrom,0.529177);
    EXPECT_EQ(ucell->tpiba,ModuleBase::TWO_PI);
    EXPECT_EQ(ucell->tpiba2,ModuleBase::TWO_PI*ModuleBase::TWO_PI);
    EXPECT_EQ(ucell->a1.x ,10.0);
    EXPECT_EQ(ucell->a2.y ,10.0);
    EXPECT_EQ(ucell->a3.z ,10.0);
    EXPECT_EQ(ucell->omega,1000.0);
    EXPECT_EQ(ucell->GT.e11,0.1);
    EXPECT_EQ(ucell->GT.e22,0.1);
    EXPECT_EQ(ucell->GT.e33,0.1);
    EXPECT_EQ(ucell->G.e11,0.1);
    EXPECT_EQ(ucell->G.e22,0.1);
    EXPECT_EQ(ucell->G.e33,0.1);

    for (int it = 0; it < ucell->ntype; it++) {
        Atom* atom = &ucell->atoms[it];
        for (int ia = 0; ia < atom->na; ia++) {
            EXPECT_EQ(atom->tau[ia].x,1);
            EXPECT_EQ(atom->tau[ia].y,1);
            EXPECT_EQ(atom->tau[ia].z,1);
        }
    }
    ofs_running.close();
    remove("setup_cell.tmp");
}


#ifdef __MPI
#include "mpi.h"
int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);

    MPI_Comm_size(MPI_COMM_WORLD,&GlobalV::NPROC);
    MPI_Comm_rank(MPI_COMM_WORLD,&GlobalV::MY_RANK);

    int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}
#endif
