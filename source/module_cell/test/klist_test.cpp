#include<iostream>
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <streambuf>
#define private public
#include "../klist.h"
#include "module_elecstate/magnetism.h"
#include "module_hamilt_pw/hamilt_pwdft/VNL_in_pw.h"
#include "module_cell/atom_pseudo.h"
#include "module_cell/atom_spec.h"
#include "module_cell/unitcell.h"
#include "module_cell/pseudo_nc.h"
#include "module_cell/setup_nonlocal.h"
#include "module_hamilt_pw/hamilt_pwdft/parallel_grid.h"
#include "module_cell/parallel_kpoints.h"
#include "module_io/berryphase.h"
#include "module_orbital/ORB_gaunt_table.h"

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
	UnitCell ucell;
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
 *     - renew and memory allocation  
 *     - read from file (Line_Cartesian kpoint file) 
 *     - setup kup and kdw after vc (different spin cases)
 */


class KlistTest : public testing::Test
{
protected:
	std::unique_ptr<K_Vectors> kv{new K_Vectors};
	std::ifstream ifs;
	std::ofstream ofs;
	std::string output;
};

TEST_F(KlistTest, Construct)
{
	EXPECT_EQ(kv->nks,0);
	EXPECT_EQ(kv->nkstot,0);
	EXPECT_EQ(kv->nkstot_ibz,0);
	EXPECT_EQ(kv->nspin,0);
	EXPECT_EQ(kv->k_nkstot,0);
	EXPECT_FALSE(kv->kc_done);
	EXPECT_FALSE(kv->kd_done);
}

TEST_F(KlistTest, MP)
{
	kv->nmp[0] = 2;
	kv->nmp[1] = 2;
	kv->nmp[2] = 2;
	kv->koffset[0] = 0;
	kv->koffset[1] = 0;
	kv->koffset[2] = 0;
	kv->nspin = 1;
	int k_type = 0;
	kv->Monkhorst_Pack(kv->nmp,kv->koffset,k_type);
	/*
	std::cout << " " <<std::endl;
	for (int ik=0;ik<kv->nkstot;ik++)
	{
		std::cout<<kv->kvec_d[ik]<<std::endl;
	}
	*/
	std::unique_ptr<K_Vectors> kv1{new K_Vectors};
	kv1->nmp[0] = 2;
	kv1->nmp[1] = 2;
	kv1->nmp[2] = 2;
	kv1->koffset[0] = 1;
	kv1->koffset[1] = 1;
	kv1->koffset[2] = 1;
	kv1->nspin = 1;
	k_type = 1;
	kv1->Monkhorst_Pack(kv1->nmp,kv1->koffset,k_type);
	//std::cout << " " <<std::endl;
	for (int ik=0;ik<kv1->nkstot;ik++)
	{
		EXPECT_EQ(kv->kvec_d[ik].x,kv1->kvec_d[ik].x);
		EXPECT_EQ(kv->kvec_d[ik].y,kv1->kvec_d[ik].y);
		EXPECT_EQ(kv->kvec_d[ik].z,kv1->kvec_d[ik].z);
		//std::cout<<kv1->kvec_d[ik]<<std::endl;
	}
}

TEST_F(KlistTest, ReadKpointsGammaOnlyLocal)
{
	GlobalV::GAMMA_ONLY_LOCAL = 1;
	std::string kfile = "KPT_GO";
	kv->nspin = 1;
	kv->read_kpoints(kfile);
	ifs.open("KPT_GO");
	std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
	EXPECT_THAT(str,testing::HasSubstr("Gamma"));
	EXPECT_THAT(str,testing::HasSubstr("1 1 1 0 0 0"));
	ifs.close();
	GlobalV::GAMMA_ONLY_LOCAL = 0; //this is important for the following tests because it is global
}

TEST_F(KlistTest, ReadKpointsKspacing)
{
	kv->nspin = 1;
	GlobalC::ucell.latvec.e11 = 10.0; GlobalC::ucell.latvec.e12 = 0.0; GlobalC::ucell.latvec.e13 = 0.0;
	GlobalC::ucell.latvec.e21 = 0.0; GlobalC::ucell.latvec.e22 = 10.0; GlobalC::ucell.latvec.e23 = 0.0;
	GlobalC::ucell.latvec.e31 = 0.0; GlobalC::ucell.latvec.e32 = 0.0; GlobalC::ucell.latvec.e33 = 10.0;
	GlobalC::ucell.GT = GlobalC::ucell.latvec.Inverse();
	GlobalC::ucell.G = GlobalC::ucell.GT.Transpose();
	GlobalC::ucell.lat0 = 1.8897261254578281;
	GlobalV::KSPACING = 0.052918; // 0.52918/Bohr = 1/A
	std::string k_file = "./support/KPT3";
	kv->read_kpoints(k_file);
	EXPECT_EQ(kv->nkstot,343);
	GlobalV::KSPACING=0.0;
}

TEST_F(KlistTest, ReadKpointsGamma)
{
	std::string k_file = "./support/KPT";
	kv->nspin = 1;
	kv->read_kpoints(k_file);
	EXPECT_EQ(kv->nkstot,512);
}

TEST_F(KlistTest, ReadKpointsMP)
{
	std::string k_file = "./support/KPT1";
	kv->nspin = 1;
	kv->read_kpoints(k_file);
	EXPECT_EQ(kv->nkstot,512);
}

TEST_F(KlistTest, ReadKpointsLine)
{
	ModuleSymmetry::Symmetry::symm_flag=0; 
	// symm_flag is required in read_kpoints for a k list
	std::string k_file = "./support/KPT2";
	kv->nspin = 1;
	kv->read_kpoints(k_file);
	EXPECT_EQ(kv->nkstot,122);
}

TEST_F(KlistTest, ReadKpointsCartesian)
{
	std::string k_file = "./support/KPT4";
        //Cartesian: non-spin case nspin=1
	kv->nspin = 1;
	kv->read_kpoints(k_file);
	EXPECT_EQ(kv->kvec_c.size(),5);
	//spin case nspin=2
	kv->nspin = 2;
	kv->read_kpoints(k_file);
	EXPECT_EQ(kv->kvec_c.size(),10);
}

TEST_F(KlistTest, ReadKpointsLineCartesian)
{
	std::string k_file = "./support/KPT5";
	//Line Cartesian: non-spin case nspin=1
	kv->nspin = 1;
	kv->set_kup_and_kdw();
	// Read from k point file under the case of Line_Cartesian.
	kv->read_kpoints(k_file);
	EXPECT_EQ(kv->nkstot,51);
	EXPECT_EQ(kv->kvec_c.size(),51);
	//Line Cartesian: spin case nspin=2
	kv->nspin = 2;
	// Read from k point file under the case of Line_Cartesian.
	kv->read_kpoints(k_file);
	EXPECT_EQ(kv->nkstot,51);
	EXPECT_EQ(kv->kvec_c.size(),102);
}

TEST_F(KlistTest, ReadKpointsDirect)
{
	std::string k_file = "./support/KPT6";
	kv->nspin = 1;
	kv->set_kup_and_kdw();
	// Read from k point file under the case of Direct
	kv->read_kpoints(k_file);
	EXPECT_EQ(kv->nkstot,6);
	EXPECT_TRUE(kv->kd_done);
}

TEST_F(KlistTest, ReadKpointsWarning1)
{
	std::string k_file = "arbitrary_1";
	kv->nspin = 1;
	GlobalV::ofs_warning.open("klist_tmp_warning_1");
	EXPECT_NO_THROW(kv->read_kpoints(k_file));
	GlobalV::ofs_warning.close();
	std::ifstream ifs;
	ifs.open("klist_tmp_warning_1");
	std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
	EXPECT_THAT(str, testing::HasSubstr("Can't find File name : arbitrary_1"));
	ifs.close();
	remove("klist_tmp_warning_1");
}

TEST_F(KlistTest, ReadKpointsWarning2)
{
	std::string k_file = "arbitrary_2";
	std::ofstream ofs;
	ofs.open(k_file.c_str());
	ofs<<"ARBITRARY";
	ofs.close();
	kv->nspin = 1;
	GlobalV::ofs_warning.open("klist_tmp_warning_2");
	EXPECT_NO_THROW(kv->read_kpoints(k_file));
	GlobalV::ofs_warning.close();
	std::ifstream ifs;
	ifs.open("klist_tmp_warning_2");
	std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
	EXPECT_THAT(str, testing::HasSubstr("symbol K_POINTS not found."));
	ifs.close();
	remove("klist_tmp_warning_2");
	remove("arbitrary_2");
}

TEST_F(KlistTest, ReadKpointsWarning3)
{
	std::string k_file = "arbitrary_3";
	std::ofstream ofs;
	ofs.open(k_file.c_str());
	ofs<<"KPOINTS"<<std::endl;
	ofs<<"100001"<<std::endl;
	ofs.close();
	kv->nspin = 1;
	GlobalV::ofs_warning.open("klist_tmp_warning_3");
	EXPECT_NO_THROW(kv->read_kpoints(k_file));
	GlobalV::ofs_warning.close();
	std::ifstream ifs;
	ifs.open("klist_tmp_warning_3");
	std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
	EXPECT_THAT(str, testing::HasSubstr("nkstot > MAX_KPOINTS"));
	ifs.close();
	remove("klist_tmp_warning_3");
	remove("arbitrary_3");
}

TEST_F(KlistTest, ReadKpointsWarning4)
{
	std::string k_file = "arbitrary_4";
	std::ofstream ofs;
	ofs.open(k_file.c_str());
	ofs<<"KPOINTS"<<std::endl;
	ofs<<"0"<<std::endl;
	ofs<<"arbitrary"<<std::endl;
	ofs.close();
	kv->nspin = 1;
	GlobalV::ofs_warning.open("klist_tmp_warning_4");
	EXPECT_NO_THROW(kv->read_kpoints(k_file));
	GlobalV::ofs_warning.close();
	std::ifstream ifs;
	ifs.open("klist_tmp_warning_4");
	std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
	EXPECT_THAT(str, testing::HasSubstr("Error: neither Gamma nor Monkhorst-Pack."));
	ifs.close();
	remove("klist_tmp_warning_4");
	remove("arbitrary_4");
}

TEST_F(KlistTest, ReadKpointsWarning5)
{
	std::string k_file = "arbitrary_5";
	std::ofstream ofs;
	ofs.open(k_file.c_str());
	ofs<<"KPOINTS"<<std::endl;
	ofs<<"100000"<<std::endl;
	ofs<<"arbitrary"<<std::endl;
	ofs.close();
        //Cartesian: non-spin case nspin=1
	kv->nspin = 1;
	GlobalV::ofs_warning.open("klist_tmp_warning_5");
	EXPECT_NO_THROW(kv->read_kpoints(k_file));
	GlobalV::ofs_warning.close();
	std::ifstream ifs;
	ifs.open("klist_tmp_warning_5");
	std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
	EXPECT_THAT(str, testing::HasSubstr("Error : neither Cartesian nor Direct kpoint"));
	ifs.close();
	remove("klist_tmp_warning_5");
	remove("arbitrary_5");
}

TEST_F(KlistTest, ReadKpointsWarning6)
{
	std::string k_file = "arbitrary_6";
	std::ofstream ofs;
	ofs.open(k_file.c_str());
	ofs<<"KPOINTS"<<std::endl;
	ofs<<"100000"<<std::endl;
	ofs<<"Line_Cartesian"<<std::endl;
	ofs.close();
        //Cartesian: non-spin case nspin=1
	kv->nspin = 1;
	ModuleSymmetry::Symmetry::symm_flag = 1;
	GlobalV::ofs_warning.open("klist_tmp_warning_6");
	EXPECT_NO_THROW(kv->read_kpoints(k_file));
	GlobalV::ofs_warning.close();
	std::ifstream ifs;
	ifs.open("klist_tmp_warning_6");
	std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
	EXPECT_THAT(str, testing::HasSubstr("Line mode of k-points is open, please set symmetry to 0 or -1"));
	ifs.close();
	remove("klist_tmp_warning_6");
	remove("arbitrary_6");
	ModuleSymmetry::Symmetry::symm_flag = 0;
}

TEST_F(KlistTest, ReadKpointsWarning7)
{
	std::string k_file = "arbitrary_7";
	std::ofstream ofs;
	ofs.open(k_file.c_str());
	ofs<<"KPOINTS"<<std::endl;
	ofs<<"100000"<<std::endl;
	ofs<<"Line_Direct"<<std::endl;
	ofs.close();
	kv->nspin = 1;
	ModuleSymmetry::Symmetry::symm_flag = 1;
	GlobalV::ofs_warning.open("klist_tmp_warning_6");
	EXPECT_NO_THROW(kv->read_kpoints(k_file));
	GlobalV::ofs_warning.close();
	std::ifstream ifs;
	ifs.open("klist_tmp_warning_6");
	std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
	EXPECT_THAT(str, testing::HasSubstr("Line mode of k-points is open, please set symmetry to 0 or -1"));
	ifs.close();
	remove("klist_tmp_warning_6");
	remove("arbitrary_6");
	ModuleSymmetry::Symmetry::symm_flag = 0;
}

TEST_F(KlistTest, SetKupKdown)
{
	std::string k_file = "./support/KPT4";
	//Cartesian: non-spin case nspin=1
	kv->nspin = 1;
	kv->read_kpoints(k_file);
	kv->set_kup_and_kdw();
	for (int ik=0;ik<5;ik++)
	{
		EXPECT_EQ(kv->isk[ik],0);

	}
	kv->nspin = 4;
	kv->read_kpoints(k_file);
	kv->set_kup_and_kdw();
	for (int ik=0;ik<5;ik++)
	{
		EXPECT_EQ(kv->isk[ik],0);
		EXPECT_EQ(kv->isk[ik+5],0);
		EXPECT_EQ(kv->isk[ik+10],0);
		EXPECT_EQ(kv->isk[ik+15],0);
	}
	kv->nspin = 2;
	kv->read_kpoints(k_file);
	kv->set_kup_and_kdw();
	for (int ik=0;ik<5;ik++)
	{
		EXPECT_EQ(kv->isk[ik],0);
		EXPECT_EQ(kv->isk[ik+5],1);
	}
}

TEST_F(KlistTest, SetKupKdownAfterVC)
{
	std::string k_file = "./support/KPT4";
	kv->nspin = 1;
	kv->read_kpoints(k_file);
	kv->set_kup_and_kdw_after_vc();
	for (int ik=0;ik<5;ik++)
	{
		EXPECT_EQ(kv->isk[ik],0);
	}
	kv->nspin = 4;
	kv->read_kpoints(k_file);
	kv->set_kup_and_kdw_after_vc();
	for (int ik=0;ik<5;ik++)
	{
		EXPECT_EQ(kv->isk[ik],0);
		EXPECT_EQ(kv->isk[ik+5],0);
		EXPECT_EQ(kv->isk[ik+10],0);
		EXPECT_EQ(kv->isk[ik+15],0);
	}
	kv->nspin = 2;
	kv->read_kpoints(k_file);
	kv->set_kup_and_kdw_after_vc();
	for (int ik=0;ik<5;ik++)
	{
		EXPECT_EQ(kv->isk[ik],0);
		EXPECT_EQ(kv->isk[ik+5],1);
	}
}

TEST_F(KlistTest, SetBothKvecAfterVC)
{
	kv->nspin = 1;
	GlobalC::ucell.latvec.e11 = 10.0; GlobalC::ucell.latvec.e12 = 0.0; GlobalC::ucell.latvec.e13 = 0.0;
	GlobalC::ucell.latvec.e21 = 0.0; GlobalC::ucell.latvec.e22 = 10.0; GlobalC::ucell.latvec.e23 = 0.0;
	GlobalC::ucell.latvec.e31 = 0.0; GlobalC::ucell.latvec.e32 = 0.0; GlobalC::ucell.latvec.e33 = 10.0;
	GlobalC::ucell.GT = GlobalC::ucell.latvec.Inverse();
	GlobalC::ucell.G = GlobalC::ucell.GT.Transpose();
	kv->nkstot = 1;
	GlobalV::ofs_running.open("tmp_klist_1");
	kv->renew(kv->nkstot);
	kv->kvec_c[0].x = 0;
	kv->kvec_c[0].y = 0;
	kv->kvec_c[0].z = 0;
	kv->set_both_kvec_after_vc(GlobalC::ucell.G,GlobalC::ucell.latvec);
	EXPECT_TRUE(kv->kd_done);
	EXPECT_TRUE(kv->kc_done);
	EXPECT_DOUBLE_EQ(kv->kvec_d[0].x,0);
	EXPECT_DOUBLE_EQ(kv->kvec_d[0].y,0);
	EXPECT_DOUBLE_EQ(kv->kvec_d[0].z,0);
	GlobalV::ofs_running.close();
	remove("tmp_klist_1");
}

TEST_F(KlistTest, PrintKlists)
{
	kv->nspin = 1;
	GlobalC::ucell.latvec.e11 = 10.0; GlobalC::ucell.latvec.e12 = 0.0; GlobalC::ucell.latvec.e13 = 0.0;
	GlobalC::ucell.latvec.e21 = 0.0; GlobalC::ucell.latvec.e22 = 10.0; GlobalC::ucell.latvec.e23 = 0.0;
	GlobalC::ucell.latvec.e31 = 0.0; GlobalC::ucell.latvec.e32 = 0.0; GlobalC::ucell.latvec.e33 = 10.0;
	GlobalC::ucell.GT = GlobalC::ucell.latvec.Inverse();
	GlobalC::ucell.G = GlobalC::ucell.GT.Transpose();
	kv->nkstot = 1;
	kv->nks = 1;
	GlobalV::ofs_running.open("tmp_klist_2");
	kv->renew(kv->nkstot);
	kv->kvec_c[0].x = 0;
	kv->kvec_c[0].y = 0;
	kv->kvec_c[0].z = 0;
	kv->set_both_kvec_after_vc(GlobalC::ucell.G,GlobalC::ucell.latvec);
	EXPECT_TRUE(kv->kd_done);
	kv->print_klists(GlobalV::ofs_running);
	GlobalV::ofs_running.close();
	remove("tmp_klist_2");
}

TEST_F(KlistTest, PrintKlistsWarnigQuit)
{
	kv->nspin = 1;
	GlobalC::ucell.latvec.e11 = 10.0; GlobalC::ucell.latvec.e12 = 0.0; GlobalC::ucell.latvec.e13 = 0.0;
	GlobalC::ucell.latvec.e21 = 0.0; GlobalC::ucell.latvec.e22 = 10.0; GlobalC::ucell.latvec.e23 = 0.0;
	GlobalC::ucell.latvec.e31 = 0.0; GlobalC::ucell.latvec.e32 = 0.0; GlobalC::ucell.latvec.e33 = 10.0;
	GlobalC::ucell.GT = GlobalC::ucell.latvec.Inverse();
	GlobalC::ucell.G = GlobalC::ucell.GT.Transpose();
	kv->nkstot = 1;
	kv->nks = 2;
	kv->renew(kv->nkstot);
	kv->kvec_c[0].x = 0;
	kv->kvec_c[0].y = 0;
	kv->kvec_c[0].z = 0;
	GlobalV::FINAL_SCF = 0;
	kv->set_both_kvec_after_vc(GlobalC::ucell.G,GlobalC::ucell.latvec);
	EXPECT_TRUE(kv->kd_done);
	testing::internal::CaptureStdout();
	EXPECT_EXIT(kv->print_klists(GlobalV::ofs_running),
			::testing::ExitedWithCode(0),"");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("nkstot < nks"));
}

TEST_F(KlistTest, SetBothKvec)
{
	kv->nspin = 1;
	GlobalC::ucell.latvec.e11 = 10.0; GlobalC::ucell.latvec.e12 = 0.0; GlobalC::ucell.latvec.e13 = 0.0;
	GlobalC::ucell.latvec.e21 = 0.0; GlobalC::ucell.latvec.e22 = 10.0; GlobalC::ucell.latvec.e23 = 0.0;
	GlobalC::ucell.latvec.e31 = 0.0; GlobalC::ucell.latvec.e32 = 0.0; GlobalC::ucell.latvec.e33 = 10.0;
	GlobalC::ucell.GT = GlobalC::ucell.latvec.Inverse();
	GlobalC::ucell.G = GlobalC::ucell.GT.Transpose();
	kv->nkstot = 1;
	kv->nks = 1;
	kv->renew(kv->nkstot);
	kv->kvec_d[0].x = 0.0;
	kv->kvec_d[0].y = 0.0;
	kv->kvec_d[0].z = 0.0;
	kv->kc_done = 0;
	kv->kd_done = 1;
	std::string skpt;
	kv->set_both_kvec(GlobalC::ucell.G,GlobalC::ucell.latvec,skpt);
	EXPECT_TRUE(kv->kc_done);
	kv->kc_done = 1;
	kv->kd_done = 0;
	kv->set_both_kvec(GlobalC::ucell.G,GlobalC::ucell.latvec,skpt);
	EXPECT_TRUE(kv->kd_done);
}

#undef private
