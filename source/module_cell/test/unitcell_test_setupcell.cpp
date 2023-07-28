#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "memory"
#include "module_base/mathzone.h"
#include "module_base/global_variable.h"
#include "module_cell/unitcell.h"
#include<vector>
#include<valarray>
#include <streambuf>
#include "prepare_unitcell.h"

#ifdef __LCAO
#include "module_basis/module_ao/ORB_read.h"
InfoNonlocal::InfoNonlocal(){}
InfoNonlocal::~InfoNonlocal(){}
LCAO_Orbitals::LCAO_Orbitals(){}
LCAO_Orbitals::~LCAO_Orbitals(){}
#endif
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

//mock function
#ifdef __LCAO
void LCAO_Orbitals::bcast_files(
	const int &ntype_in,
	const int &my_rank)
{
	return;
}

class UcellTest : public ::testing::Test
{
protected:
	std::unique_ptr<UnitCell> ucell{new UnitCell};
	std::string output;
	void SetUp()
{
    	ucell->lmaxmax = 2;
}
};

using UcellDeathTest = UcellTest;

TEST_F(UcellTest,SetupCellS1)
{
	std::string fn = "./support/STRU_MgO";
	std::ofstream ofs_running;
	ofs_running.open("setup_cell.tmp");
	GlobalV::NSPIN = 1;
	ucell->ntype = 2;
	ucell->setup_cell(fn,ofs_running);
	ofs_running.close();
	remove("setup_cell.tmp");
}

TEST_F(UcellTest,SetupCellS2)
{
	std::string fn = "./support/STRU_MgO";
	std::ofstream ofs_running;
	ofs_running.open("setup_cell.tmp");
	GlobalV::NSPIN = 2;
	ucell->ntype = 2;
	ucell->setup_cell(fn,ofs_running);
	ofs_running.close();
	remove("setup_cell.tmp");
}

TEST_F(UcellTest,SetupCellS4)
{
	std::string fn = "./support/STRU_MgO";
	std::ofstream ofs_running;
	ofs_running.open("setup_cell.tmp");
	GlobalV::NSPIN = 4;
	ucell->ntype = 2;
	ucell->setup_cell(fn,ofs_running);
	ofs_running.close();
	remove("setup_cell.tmp");
}

TEST_F(UcellDeathTest,SetupCellWarning1)
{
	std::string fn = "./STRU_MgO";
	std::ofstream ofs_running;
	ofs_running.open("setup_cell.tmp");
	ucell->ntype = 2;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(ucell->setup_cell(fn,ofs_running),::testing::ExitedWithCode(0),"");
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
	ucell->ntype = 2;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(ucell->setup_cell(fn,ofs_running),::testing::ExitedWithCode(0),"");
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
	GlobalV::NSPIN = 1;
	ucell->ntype = 2;
	delete[] ucell->magnet.start_magnetization;
	ucell->magnet.start_magnetization = new double[ucell->ntype];
	ucell->setup_cell(fn,ofs_running);
	ucell->setup_cell_after_vc(ofs_running);
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
#endif
