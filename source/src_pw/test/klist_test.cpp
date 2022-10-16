#include "gtest/gtest.h"
#include "src_pw/magnetism.h"
#include "src_pw/VNL_in_pw.h"
#include "module_cell/atom_pseudo.h"
#include "module_cell/atom_spec.h"
#include "module_cell/unitcell_pseudo.h"
#include "module_cell/pseudo_nc.h"
#include "module_cell/setup_nonlocal.h"
#include "src_parallel/parallel_grid.h"
#include "src_parallel/parallel_kpoints.h"
#include "src_io/berryphase.h"
#include "module_orbital/ORB_gaunt_table.h"
#include<iostream>

bool berryphase::berry_phase_flag=0;

pseudo_nc::pseudo_nc(){}
pseudo_nc::~pseudo_nc(){}
Atom::Atom(){}
Atom::~Atom(){}
Atom_pseudo::Atom_pseudo(){}
Atom_pseudo::~Atom_pseudo(){}
InfoNonlocal::InfoNonlocal(){}
InfoNonlocal::~InfoNonlocal(){}
UnitCell::UnitCell(){}
UnitCell::~UnitCell(){}
UnitCell_pseudo::UnitCell_pseudo(){}
UnitCell_pseudo::~UnitCell_pseudo(){}
Magnetism::Magnetism(){}
Magnetism::~Magnetism(){}
ORB_gaunt_table::ORB_gaunt_table(){}
ORB_gaunt_table::~ORB_gaunt_table(){}
pseudopot_cell_vl::pseudopot_cell_vl(){}
pseudopot_cell_vl::~pseudopot_cell_vl(){}
pseudopot_cell_vnl::pseudopot_cell_vnl(){}
pseudopot_cell_vnl::~pseudopot_cell_vnl(){}

namespace GlobalC
{
	Parallel_Kpoints Pkpoints;
	UnitCell_pseudo ucell;
}

/************************************************
 *  unit test of class K_Vectors
 ***********************************************/

/**
 * - Tested Functions:
 *   - K_Vectors()
 *     - basic parameters (nks,nkstot,nkstot_ibz) are set
 *   - read_kpoints()
 *     - read from file
 *     - generate KPT from kspacing parameter
 */

#define private public
#include "src_pw/klist.h"

class KlistTest : public testing::Test
{
protected:
	K_Vectors kv;
};

TEST_F(KlistTest, Construct)
{
	EXPECT_EQ(kv.nks,0);
	EXPECT_EQ(kv.nkstot,0);
	EXPECT_EQ(kv.nkstot_ibz,0);
	EXPECT_EQ(kv.nspin,0);
	EXPECT_EQ(kv.k_nkstot,0);
	EXPECT_FALSE(kv.kc_done);
	EXPECT_FALSE(kv.kd_done);
}

TEST_F(KlistTest, MP)
{
	kv.nmp[0] = 2;
	kv.nmp[1] = 2;
	kv.nmp[2] = 2;
	kv.koffset[0] = 0;
	kv.koffset[1] = 0;
	kv.koffset[2] = 0;
	kv.nspin = 1;
	int k_type = 0;
	kv.Monkhorst_Pack(kv.nmp,kv.koffset,k_type);
	/*
	std::cout << " " <<std::endl;
	for (int ik=0;ik<kv.nkstot;ik++)
	{
		std::cout<<kv.kvec_d[ik]<<std::endl;
	}
	*/
	K_Vectors kv1;
	kv1.nmp[0] = 2;
	kv1.nmp[1] = 2;
	kv1.nmp[2] = 2;
	kv1.koffset[0] = 1;
	kv1.koffset[1] = 1;
	kv1.koffset[2] = 1;
	kv1.nspin = 1;
	k_type = 1;
	kv1.Monkhorst_Pack(kv1.nmp,kv1.koffset,k_type);
	//std::cout << " " <<std::endl;
	for (int ik=0;ik<kv1.nkstot;ik++)
	{
		EXPECT_EQ(kv.kvec_d[ik].x,kv1.kvec_d[ik].x);
		EXPECT_EQ(kv.kvec_d[ik].y,kv1.kvec_d[ik].y);
		EXPECT_EQ(kv.kvec_d[ik].z,kv1.kvec_d[ik].z);
		//std::cout<<kv1.kvec_d[ik]<<std::endl;
	}
}

TEST_F(KlistTest, ReadMesh)
{
	std::string k_file = "KPT";
	kv.nspin = 1;
	kv.read_kpoints(k_file);
	EXPECT_EQ(kv.nkstot,512);
	remove("KPT");
}

TEST_F(KlistTest, ReadMP)
{
	K_Vectors kv;
	std::string k_file = "KPT1";
	kv.nspin = 1;
	kv.read_kpoints(k_file);
	EXPECT_EQ(kv.nkstot,512);
}

TEST_F(KlistTest, ReadList)
{
	ModuleSymmetry::Symmetry::symm_flag=0; 
	// symm_flag is required in read_kpoints for a k list
	K_Vectors kv;
	std::string k_file = "KPT2";
	kv.nspin = 1;
	kv.read_kpoints(k_file);
	EXPECT_EQ(kv.nkstot,122);
	remove("KPT2");
}

TEST_F(KlistTest, Kspacing)
{
	kv.nspin = 1;
	GlobalC::ucell.latvec.e11 = 10.0; GlobalC::ucell.latvec.e12 = 0.0; GlobalC::ucell.latvec.e13 = 0.0;
	GlobalC::ucell.latvec.e21 = 0.0; GlobalC::ucell.latvec.e22 = 10.0; GlobalC::ucell.latvec.e23 = 0.0;
	GlobalC::ucell.latvec.e31 = 0.0; GlobalC::ucell.latvec.e32 = 0.0; GlobalC::ucell.latvec.e33 = 10.0;
	GlobalC::ucell.GT = GlobalC::ucell.latvec.Inverse();
	GlobalC::ucell.G = GlobalC::ucell.GT.Transpose();
	GlobalC::ucell.lat0 = 1.8897261254578281;
	GlobalV::KSPACING = 0.052918; // 0.52918/Bohr = 1/A
	std::string k_file = "KPT3";
	kv.read_kpoints(k_file);
	EXPECT_EQ(kv.nkstot,343);
	remove("KPT3");
}

#undef private
