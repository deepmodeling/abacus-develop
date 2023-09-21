#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include<streambuf>
#ifdef __MPI
#include "mpi.h"
#endif

/************************************************
 *  unit test of Atom_pseudo
 ***********************************************/

/**
 * - Tested Functions:
 *   - Atom()
 *     - constructor of class Atom
 *   - ~Atom()
 *     - deconstructor of class Atom
 *   - print_Atom()
 *     - print atomic info to file
 *   - set_index()
 *     - index between numerical atomic orbtials and quantum numbers L,m
 *   - bcast_atom()
 *     - bcast basic atomic info to all processes
 *   - bcast_atom2()
 *     - bcast norm-conserving pseudopotential info to all processes
 */

#define private public
#include "module_cell/read_pp.h"
#include "module_cell/pseudo.h"
#include "module_cell/atom_pseudo.h"
#include "module_cell/atom_spec.h"

class AtomSpecTest : public testing::Test
{
protected:
	Atom atom;
	Pseudopot_upf upf;
	std::ofstream ofs;
	std::ifstream ifs;
};

TEST_F(AtomSpecTest, PrintAtom)
{
#ifdef __MPI
	if(GlobalV::MY_RANK==0)
	{
#endif
	ofs.open("tmp_atom_info");
	atom.label = "C";
	atom.type = 1;
	atom.na = 2;
	atom.nwl = 2;
	atom.Rcut = 1.1;
	atom.nw = 14;
	atom.stapos_wf = 0;
	atom.mass = 12.0;
	delete[] atom.tau;
	atom.tau = new ModuleBase::Vector3<double>[atom.na];
	atom.tau[0].x = 0.2;
	atom.tau[0].y = 0.2;
	atom.tau[0].z = 0.2;
	atom.tau[1].x = 0.4;
	atom.tau[1].y = 0.4;
	atom.tau[1].z = 0.4;
	atom.print_Atom(ofs);
	ofs.close();
	ifs.open("tmp_atom_info");
	std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
	   	EXPECT_THAT(str, testing::HasSubstr("label = C"));
    	EXPECT_THAT(str, testing::HasSubstr("type = 1"));
    	EXPECT_THAT(str, testing::HasSubstr("na = 2"));
    	EXPECT_THAT(str, testing::HasSubstr("nwl = 2"));
    	EXPECT_THAT(str, testing::HasSubstr("Rcut = 1.1"));
    	EXPECT_THAT(str, testing::HasSubstr("nw = 14"));
    	EXPECT_THAT(str, testing::HasSubstr("stapos_wf = 0"));
    	EXPECT_THAT(str, testing::HasSubstr("mass = 12"));
    	EXPECT_THAT(str, testing::HasSubstr("atom_position(cartesian) Dimension = 2"));
	ifs.close();
	remove("tmp_atom_info");
#ifdef __MPI
	}
#endif
}

TEST_F(AtomSpecTest, SetIndex)
{
#ifdef __MPI
	if(GlobalV::MY_RANK==0)
	{
#endif
	atom.nw = 0;
	atom.nwl = 1;
	delete[] atom.l_nchi;
	atom.l_nchi = new int[atom.nwl+1];
	atom.l_nchi[0] = 2; // l:0, N:2 (arbitrary)
	atom.nw += 1*atom.l_nchi[0]; // m = 2*0+1 = 1
	atom.l_nchi[1] = 4; // l:1, N:4 (arbitrary)
	atom.nw += 3*atom.l_nchi[1]; // m = 2*1+1 = 3
	atom.set_index();
	EXPECT_EQ(atom.iw2l[13],1);
	EXPECT_EQ(atom.iw2n[13],3);
	EXPECT_EQ(atom.iw2m[13],2);
	EXPECT_EQ(atom.iw2_ylm[13],3);
	EXPECT_TRUE(atom.iw2_new[11]);
	// here is the table:
	// nw = 2 + 3*4 = 14
	//    L N m  L*L+m
	// 0  0 0 0    0
	// 1  0 1 0    0
	// 2  1 0 0    1
	// 3  1 0 1    2
	// 4  1 0 2    3
	// 5  1 1 0    1
	// 6  1 1 1    2
	// 7  1 1 2    3
	// 8  1 2 0    1
	// 9  1 2 1    2
	// 10 1 2 2    3
	// 11 1 3 0    1
	// 12 1 3 1    2
	// 13 1 3 2    3
#ifdef __MPI
	}
#endif
}

#ifdef __MPI
TEST_F(AtomSpecTest, BcastAtom)
{
	GlobalV::test_atom = 1;
	if(GlobalV::MY_RANK==0)
	{
		atom.label = "C";
		atom.type = 1;
		atom.na = 2;
		atom.nw = 0;
		atom.nwl = 1;
		atom.Rcut = 1.1;
		delete[] atom.l_nchi;
		atom.l_nchi = new int[atom.nwl+1];
		atom.l_nchi[0] = 2;
		atom.nw += atom.l_nchi[0];
		atom.l_nchi[1] = 4;
		atom.nw += 3*atom.l_nchi[1];
		atom.stapos_wf = 0;
		atom.mass = 12.0;
		delete[] atom.tau;
		delete[] atom.taud;
		delete[] atom.vel;
		delete[] atom.mag;
		delete[] atom.angle1;
		delete[] atom.angle2;
		delete[] atom.m_loc_;
		delete[] atom.mbl;
		atom.tau = new ModuleBase::Vector3<double>[atom.na];
		atom.taud = new ModuleBase::Vector3<double>[atom.na];
		atom.vel = new ModuleBase::Vector3<double>[atom.na];
		atom.mag = new double[atom.na];
		atom.angle1 = new double[atom.na];
		atom.angle2 = new double[atom.na];
		atom.m_loc_ = new ModuleBase::Vector3<double>[atom.na];
		atom.mbl = new ModuleBase::Vector3<int>[atom.na];
		atom.tau[0].x = 0.2;
		atom.tau[0].y = 0.2;
		atom.tau[0].z = 0.2;
		atom.tau[1].x = 0.4;
		atom.tau[1].y = 0.4;
		atom.tau[1].z = 0.4;
	}
	atom.bcast_atom();
	if(GlobalV::MY_RANK!=0)
	{
		EXPECT_EQ(atom.label,"C");
		EXPECT_EQ(atom.type,1);
		EXPECT_EQ(atom.na,2);
		EXPECT_EQ(atom.nwl,1);
		EXPECT_DOUBLE_EQ(atom.Rcut,1.1);
		EXPECT_EQ(atom.nw,14);
		EXPECT_EQ(atom.stapos_wf,0);
		EXPECT_DOUBLE_EQ(atom.mass,12.0);
		EXPECT_DOUBLE_EQ(atom.tau[0].x,0.2);
		EXPECT_DOUBLE_EQ(atom.tau[1].z,0.4);
	}
}

TEST_F(AtomSpecTest, BcastAtom2)
{
	if(GlobalV::MY_RANK==0)
	{
		ifs.open("./support/C.upf");
		GlobalV::PSEUDORCUT = 15.0;
		upf.read_pseudo_upf201(ifs);
		atom.ncpp.set_pseudo_nc(upf);
		ifs.close();
		EXPECT_TRUE(atom.ncpp.has_so);
	}
	atom.bcast_atom2();
	if(GlobalV::MY_RANK!=0)
	{
		EXPECT_EQ(atom.ncpp.nbeta,6);
		EXPECT_EQ(atom.ncpp.nchi,3);
		EXPECT_DOUBLE_EQ(atom.ncpp.rho_atc[0],8.7234550809E-01);
	}
}

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
#undef private
