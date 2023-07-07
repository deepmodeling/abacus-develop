#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "./for_test_sltk_atom_arrange.h"
#include "../sltk_atom_arrange.h"

#include<iostream>
#include<string>

/************************************************
 *  unit test of atom_arrange
 ***********************************************/

/**
 * - Tested Functions:
 *   - atom_arrange::delete_vector(void)
 *     - delete vector
 *   - atom_arrange::set_sr_NL
 * 	   - set the sr: search radius including nonlocal beta
 */

class SltkAtomArrangeTest : public testing::Test
{
};

TEST_F(SltkAtomArrangeTest, setsrNL)
{
    atom_arrange test;
    std::string teststring="m";
    double rcutmax_Phi=1; 
    double rcutmax_Beta=2; 
    bool gamma_only_local=true;
    double test_sr=0;
    std::ofstream ofs;
    ofs.open("./to_test_arrange.txt");
    test_sr=test.set_sr_NL(ofs,teststring,rcutmax_Phi,rcutmax_Beta,gamma_only_local);
    EXPECT_DOUBLE_EQ(test_sr,2.01);

    gamma_only_local=false;
    test_sr=test.set_sr_NL(ofs,teststring,rcutmax_Phi,rcutmax_Beta,gamma_only_local);
    EXPECT_DOUBLE_EQ(test_sr,6.01);

    teststring="no";
    test_sr=test.set_sr_NL(ofs,teststring,rcutmax_Phi,rcutmax_Beta,gamma_only_local);
    ofs.close();
    std::ifstream ifs;
    std::string test2,s;
    ifs.open("./to_test_arrange.txt");
    while (getline(ifs,s))
    {
        test2 += s;
    }
    std::string test3=" >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
	" |                                                                    |"
	" | Search adjacent atoms:                                             |"
	" | Set the adjacent atoms for each atom and set the periodic boundary |"
	" | condition for the atoms on real space FFT grid. For k-dependent    |"
	" | algorithm, we also need to set the sparse H and S matrix element   |"
	" | for each atom.                                                     |"
	" |                                                                    |"
	" <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
    " SETUP SEARCHING RADIUS FOR PROGRAM TO SEARCH ADJACENT ATOMS"
    "                  longest orb rcut (Bohr) = 1   longest nonlocal projector rcut (Bohr) = 2";
    EXPECT_EQ(test2,test3);
	
    ifs.close();
    remove("./to_test_arrange");
}

TEST_F(SltkAtomArrangeTest, DeleteVector)
{
	atom_arrange test;
	std::ofstream ofs("./testforarrange.txt");
	bool pbc_flag = true;
	Grid_Driver grid_d(1,2,3); 
	UnitCell ucell; 
	double search_radius_bohr=2; 
	int test_atom_in= 0;
	grid_d.init_cell_flag=false;

	int dx=1,dy=1,dz=1;
	grid_d.Cell = new CellSet**[dx];
	for (int i = 0;i < dx;i++)
	{
		grid_d.Cell[i] = new CellSet*[dy];
		
		for (int j = 0;j < dy;j++)
		{
			grid_d.Cell[i][j] = new CellSet[dz];
			
			for (int k = 0;k < dz;k++)
			{
				grid_d.Cell[i][j][k].length = 0;
			}
		}
	}
	std::string output;
	testing::internal::CaptureStdout();
	test.delete_vector(ofs,pbc_flag,grid_d,ucell,search_radius_bohr,test_atom_in);
	output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output,testing::HasSubstr("test"));
	ofs.close();
	
}

TEST_F(SltkAtomArrangeTest, DeleteVector2)
{
	std::ofstream ofs("./testforarrange.txt");
	atom_arrange test2;
	bool pbc_flag = true;
	Grid_Driver grid_d2(1,2,3); 
	UnitCell ucell; 
	double search_radius_bohr=2; 
	int test_atom_in= 0;
	
	int dx=1,dy=1,dz=1;
	grid_d2.dx=1;
	grid_d2.dy=1;
	grid_d2.dz=1;
	grid_d2.init_cell_flag=true;
	grid_d2.Cell = new CellSet**[grid_d2.dx];
	for (int i = 0;i < grid_d2.dx;i++)
	{
		grid_d2.Cell[i] = new CellSet*[grid_d2.dy];
		
		for (int j = 0;j < grid_d2.dy;j++)
		{
			grid_d2.Cell[i][j] = new CellSet[grid_d2.dz];
			
			for (int k = 0;k < grid_d2.dz;k++)
			{
				grid_d2.Cell[i][j][k].length = 0;
			}
		}
	}
	EXPECT_EQ(grid_d2.Cell[0][0][0].length,0);
	std::string output;
	testing::internal::CaptureStdout();
	test2.delete_vector(ofs,pbc_flag,grid_d2,ucell,search_radius_bohr,test_atom_in);
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("test"));
//	EXPECT_EQ(grid_d.Cell,NULL);
	EXPECT_FALSE(grid_d2.init_cell_flag);
	ofs.close();
}


