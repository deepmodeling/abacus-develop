#include "gtest/gtest.h"
#include "gmock/gmock.h"
#define private public
#include "module_parameter/parameter.h"
#undef private
#include "memory"
#include "module_base/mathzone.h"
#include "module_base/global_variable.h"
#include "module_cell/unitcell.h"
#include "module_cell/read_stru.h"
#include<vector>
#include<valarray>

Magnetism::Magnetism()
{
	this->tot_magnetization = 0.0;
	this->abs_magnetization = 0.0;
	this->start_magnetization = nullptr;
}
Magnetism::~Magnetism()
{
	delete[] this->start_magnetization;
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
	void SetUp()
    {
    	ucell->lmaxmax = 2;
        ucell->ntype   = 2;
        ucell->atom_mass.resize(ucell->ntype);
        ucell->atom_label.resize(ucell->ntype);
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
	PARAM.input.test_pseudo_cell = 2;
	EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running,*ucell));
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
	PARAM.input.test_pseudo_cell = 2;
	PARAM.input.basis_type = "pw";
	//call read_atom_species
	EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running,*ucell));
	EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat));
	EXPECT_DOUBLE_EQ(ucell->latvec.e11,4.27957);
	EXPECT_DOUBLE_EQ(ucell->latvec.e22,4.27957);
	EXPECT_DOUBLE_EQ(ucell->latvec.e33,4.27957);
	//mandatory preliminaries
	delete[] ucell->magnet.start_magnetization;
	ucell->magnet.start_magnetization = new double[ucell->ntype];
	//call read_atom_positions
	EXPECT_NO_THROW(ucell->read_atom_positions(ifa,ofs_running,ofs_warning));
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
	PARAM.input.nspin = 1;
	ucell->setup_cell(fn,ofs_running);
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
