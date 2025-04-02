#include "gtest/gtest.h"
#include "gmock/gmock.h"
#define private public
#include "module_parameter/parameter.h"
#undef private
#include "memory"
#include "module_base/mathzone.h"
#include "module_base/global_variable.h"
#include "module_cell/unitcell.h"
#include<vector>
#include<valarray>
#include <streambuf>
#include "prepare_unitcell.h"
#include "module_cell/update_cell.h"
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
		ucell->ntype   = 2;
        ucell->atom_mass.resize(ucell->ntype);
        ucell->atom_label.resize(ucell->ntype);
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
	PARAM.input.nspin = 1;
	
	ucell->setup_cell(fn,ofs_running);
	ofs_running.close();
	remove("setup_cell.tmp");
}

TEST_F(UcellTest,SetupCellS2)
{
	std::string fn = "./support/STRU_MgO";
	std::ofstream ofs_running;
	ofs_running.open("setup_cell.tmp");
	PARAM.input.nspin = 2;
	
	ucell->setup_cell(fn,ofs_running);
	ofs_running.close();
	remove("setup_cell.tmp");
}

TEST_F(UcellTest,SetupCellS4)
{
	std::string fn = "./support/STRU_MgO";
	std::ofstream ofs_running;
	ofs_running.open("setup_cell.tmp");
	PARAM.input.nspin = 4;
	
	ucell->setup_cell(fn,ofs_running);
	ofs_running.close();
	remove("setup_cell.tmp");
}

TEST_F(UcellDeathTest,SetupCellWarning1)
{
	std::string fn = "./STRU_MgO";
	std::ofstream ofs_running;
	ofs_running.open("setup_cell.tmp");
	
	testing::internal::CaptureStdout();
	EXPECT_EXIT(ucell->setup_cell(fn,ofs_running),::testing::ExitedWithCode(1),"");
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
	EXPECT_EXIT(ucell->setup_cell(fn,ofs_running),::testing::ExitedWithCode(1),"");
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
	PARAM.input.nspin = 1;
	
	delete[] ucell->magnet.start_magnetization;
	ucell->magnet.start_magnetization = new double[ucell->ntype];

	
	ucell->setup_cell(fn,ofs_running);
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
	
	unitcell::setup_cell_after_vc(*ucell,ofs_running);
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
#endif
