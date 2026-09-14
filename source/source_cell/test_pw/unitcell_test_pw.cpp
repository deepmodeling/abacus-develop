#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "memory"
#include "source_base/mathzone.h"
#include "source_base/global_variable.h"
#include "source_cell/unitcell.h"
#include "source_cell/read_stru.h"
#include<vector>
#include<valarray>

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
 *   - ReadAtomSpecies
 *     - read_atom_species(): read header part of orbital file
 *   - ReadAtomPositions
 *     - read_atom_positions(): read atomic coordinates, velocities, magmoms
 *   - SetupCell
 *     - setup_cell(): the pw version
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

TEST_F(UcellTest,ReadAtomSpecies)
{
#ifdef __MPI
if(GlobalV::MY_RANK==0)
{
#endif
    std::string fn = "./support/STRU_MgO";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    ofs_running.open("read_atom_species.tmp");
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc, onsite_radius, deepks_setorb, rpa));
    EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat));
    EXPECT_DOUBLE_EQ(ucell->latvec.e11,4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e22,4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e33,4.27957);
    ofs_running.close();
    ifa.close();
    remove("read_atom_species.tmp");
#ifdef __MPI
}
#endif
}

TEST_F(UcellTest,ReadAtomPositions)
{
#ifdef __MPI
if(GlobalV::MY_RANK==0)
{
#endif
    std::string fn = "./support/STRU_MgO";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    std::ofstream ofs_warning;
    ofs_running.open("read_atom_species.tmp");
    ofs_warning.open("read_atom_species.warn");
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    const int nspin = 1;
    //call read_atom_species
    EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc, onsite_radius, deepks_setorb, rpa));
    EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat));
    EXPECT_DOUBLE_EQ(ucell->latvec.e11,4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e22,4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e33,4.27957);
    //call read_atom_positions
    EXPECT_NO_THROW(unitcell::read_atom_positions(*ucell, ifa, ofs_running, ofs_warning, nspin,
        basis_type, orbital_dir, init_wfc, onsite_radius, fixed_atoms, noncolin,
        calculation, esolver_type, 0));
    ofs_running.close();
    ofs_warning.close();
    ifa.close();
    remove("read_atom_species.tmp");
    remove("read_atom_species.warn");
#ifdef __MPI
}
#endif
}

TEST_F(UcellTest,SetupCell)
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
