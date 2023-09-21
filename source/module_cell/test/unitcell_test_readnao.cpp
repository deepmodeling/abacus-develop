#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "memory"
#include "module_base/mathzone.h"
#include "module_base/global_variable.h"
#include "module_cell/unitcell.h"
#include <vector>
#include <valarray>
#ifdef __MPI
#include "mpi.h"
#endif
#include "prepare_unitcell.h"

#ifdef __LCAO
InfoNonlocal::InfoNonlocal(){}
InfoNonlocal::~InfoNonlocal(){}
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
 *   - ReadNaoFlzTest
 *     - read_nao_flz(): read numerical orbital file and store in atom member variables rgrid, flz and n_rgrid
 * - NOTE: This function is only a templorary choice. Once module_nao is decoupled from old class ORB, this function will be removed.
 *         Currently the function to test is for reading-in detailed numerical orbital information for so-called lcao_in_pw calculation
 *         but now which is actually, the initialization by numerical orbitals on psi, in planewave representation
 * - To remove:
 *   - remove definition and declaration in read_atoms.cpp and unitcell.h
 *   - remove present test in CMakeLists.txt
 *   - delete this file
 */

//mock function
#ifdef __LCAO
void LCAO_Orbitals::bcast_files(
	const int &ntype_in,
	const int &my_rank)
{
	return;
}
#endif

class UcellTest : public ::testing::Test
{
protected:
    /* Here mock one unitcell object */

    UnitCell ucell;
	std::ofstream ofs;

    int error = 0;
	void SetUp()
	{
		ofs.open("running.log");

		GlobalV::global_out_dir = "./";
		GlobalV::LSPINORB = false;

		GlobalV::PSEUDORCUT = 15.0;
		GlobalV::DFT_FUNCTIONAL = "default";
		GlobalV::test_unitcell = 1;
		GlobalV::test_pseudo_cell = 1;
		GlobalV::NSPIN = 1;
		GlobalV::BASIS_TYPE = "pw";
        
        ucell.atoms = new Atom[1];
        ucell.atoms[0].label = "C";
        ucell.atoms[0].na = 1;
        ucell.atom_label = new std::string[1];
        ucell.atom_label[0] = "C";
        this->error = std::system("cp ../../../../source/module_cell/test/support/C_gga_8au_100Ry_2s2p1d.orb ./C_gga_8au_100Ry_2s2p1d.orb");
        this->error = std::system("cp ../../../../source/module_cell/test/support/C_gga_8au_100Ry_2s2p1d.orb ./H_gga_8au_100Ry_2s1p.orb");
	}
	void TearDown()
	{
		ofs.close();
        this->error = std::system("rm ./C_gga_8au_100Ry_2s2p1d.orb");
        this->error = std::system("rm ./H_gga_8au_100Ry_2s1p.orb");
	}
};

TEST_F(UcellTest, ReadNaoFlzTest)
{
    std::string fn = "C_gga_8au_100Ry_2s2p1d.orb";
    // the ichi goes from 1st s, 2nd s, 1st p, 2nd p, 1st d
    ucell.read_orb_file(0, fn, ofs, &(ucell.atoms[0]));
    ucell.read_nao_flz(0, fn, ofs, &(ucell.atoms[0]));

    EXPECT_EQ(ucell.atoms[0].n_rgrid[0], 801);
    EXPECT_DOUBLE_EQ(ucell.atoms[0].rgrid[0][0], 0);
    EXPECT_DOUBLE_EQ(ucell.atoms[0].rgrid[0][1], 0.01);
    EXPECT_DOUBLE_EQ(ucell.atoms[0].rgrid[0][800], 8);

    EXPECT_EQ(ucell.atoms[0].n_rgrid[1], 801);
    EXPECT_DOUBLE_EQ(ucell.atoms[0].rgrid[1][0], 0);
    EXPECT_DOUBLE_EQ(ucell.atoms[0].rgrid[1][1], 0.01);
    EXPECT_DOUBLE_EQ(ucell.atoms[0].rgrid[1][800], 8);

    // data of 1st s
    EXPECT_DOUBLE_EQ(ucell.atoms[0].flz[0][0], 5.368426038998e-01);
    EXPECT_DOUBLE_EQ(ucell.atoms[0].flz[0][4], 5.386156949024e-01);
    EXPECT_DOUBLE_EQ(ucell.atoms[0].flz[0][800], 0);

    // data of 2nd s
    EXPECT_DOUBLE_EQ(ucell.atoms[0].flz[1][0], -6.134205291735e-02);
    EXPECT_DOUBLE_EQ(ucell.atoms[0].flz[1][4], -5.862821438615e-02);
    EXPECT_DOUBLE_EQ(ucell.atoms[0].flz[1][800], 0);
}


int main(int argc, char **argv)
{
	testing::InitGoogleTest(&argc, argv);

	int result = RUN_ALL_TESTS();
	return result;
}
