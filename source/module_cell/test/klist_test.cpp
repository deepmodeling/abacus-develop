#include<iostream>
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <streambuf>
#include "module_base/mathzone.h"
#include "module_cell/parallel_kpoints.h"
#include "module_base/parallel_global.h"
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
#include "module_basis/module_ao/ORB_gaunt_table.h"

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
 *     - ReadKpointsGammaOnlyLocal: GlobalV::GAMMA_ONLY_LOCAL = 1
 *     - ReadKpointsKspacing: generate KPT from kspacing parameter 
 *     - ReadKpointsGamma: "Gamma" mode of `KPT` file
 *     - ReadKpointsMP: "MP" mode of `KPT` file
 *     - ReadKpointsLine: "Line" mode of `KPT` file
 *     - ReadKpointsCartesian: "Cartesian" mode of `KPT` file
 *     - ReadKpointsLineCartesian: "Line_Cartesian" mode of `KPT` file
 *     - ReadKpointsDirect: "Direct" mode of `KPT` file
 *     - ReadKpointsWarning1: Can't find KPT file through string k_file
 *     - ReadKpointsWarning2: symbol K_POINTS not found at the 1st line
 *     - ReadKpointsWarning3: nkstot > MAX_KPOINTS (100000)
 *     - ReadKpointsWarning4: Error: neither Gamma nor Monkhorst-Pack mode when nkstot=0
 *     - ReadKpointsWarning5: Error : neither Cartesian nor Direct kpoint in line mode(nkstot>0)
 *     - ReadKpointsWarning6: Line mode (Line_Cartesian) of k-points is open, please set symmetry to 0 or -1
 *     - ReadKpointsWarning7: Line mode (Line_Direct) of k-points is open, please set symmetry to 0 or -1
 *   - set_kup_and_kdw()
 *     - SetKupKdown: set basic kpoints info: kvec_c, kvec_d, wk, isk, nks, nkstot
 *       according to different spin case
 *   - set_kup_and_kdw_after_vc()
 *     - SetKupKdownAfterVC: set basic kpoints info: kvec_c, kvec_d, wk, isk, nks, nkstot
 *       according to different spin case after variable-cell optimization
 *   - set_both_kvec()
 *     - SetBothKvec: set kvec_c (cartesian coor.) and kvec_d (direct coor.)
 *     - SetBothKvecFinalSCF: same as above, with GlobalV::FINAL_SCF=1
 *   - set_both_kvec_after_vc()
 *     - SetBothKvecAfterVC: set kvec_c (cartesian coor.) and kvec_d (direct coor.)
 *       after variable-cell relaxation
 *   - print_klists()
 *     - PrintKlists: print kpoints coordinates
 *     - PrintKlistsWarningQuit: for nkstot < nks error
 *   - normalize_wk()
 *     - NormalizeWk: normalize weight of kpoints
 *   - update_use_ibz()
 *     - UpdateUseIbz: update kpoints info by ibz kpoints info
 *   - ibz_kpoint()
 *     - IbzKpoint: generate IBZ kpoints
 *     - IbzKpointIsMP: generate IBZ kpoints for non-symmetry 
 *       and Monkhorst-Pack case
 */

//abbriviated from module_symmetry/test/symmetry_test.cpp
struct atomtype_
{
    std::string atomname;
    std::vector<std::vector<double>> coordinate;
};

struct stru_
{
    int ibrav;
    std::string point_group; // Schoenflies symbol
    std::string point_group_hm; // Hermann-Mauguin notation.
    std::string space_group;
    std::vector<double> cell;
    std::vector<atomtype_> all_type;
};

std::vector<stru_> stru_lib{
    stru_{1,
          "O_h",
          "m-3m",
          "Pm-3m",
          std::vector<double>{1., 0., 0., 0., 1., 0., 0., 0., 1.},
          std::vector<atomtype_>{atomtype_{"C",
                                           std::vector<std::vector<double>>{
                                               {0., 0., 0.},
                                           }}}}
};
//used to construct cell and analyse its symmetry


class KlistTest : public testing::Test
{
protected:
	std::unique_ptr<K_Vectors> kv{new K_Vectors};
	std::ifstream ifs;
	std::ofstream ofs;
	std::ofstream ofs_running;
	std::string output;

	//used to construct cell and analyse its symmetry
	UnitCell ucell;
	void construct_ucell(stru_ &stru)
	{
	    std::vector<atomtype_> coord = stru.all_type;
	    ucell.a1 = ModuleBase::Vector3<double>(stru.cell[0], stru.cell[1], stru.cell[2]);
	    ucell.a2 = ModuleBase::Vector3<double>(stru.cell[3], stru.cell[4], stru.cell[5]);
	    ucell.a3 = ModuleBase::Vector3<double>(stru.cell[6], stru.cell[7], stru.cell[8]);
	    ucell.ntype = stru.all_type.size();
	    ucell.atoms = new Atom[ucell.ntype];
	    ucell.nat = 0;
	    ucell.latvec.e11 = ucell.a1.x; ucell.latvec.e12 = ucell.a1.y; ucell.latvec.e13 = ucell.a1.z;
	    ucell.latvec.e21 = ucell.a2.x; ucell.latvec.e22 = ucell.a2.y; ucell.latvec.e23 = ucell.a2.z;
	    ucell.latvec.e31 = ucell.a3.x; ucell.latvec.e32 = ucell.a3.y; ucell.latvec.e33 = ucell.a3.z;
	    ucell.GT = ucell.latvec.Inverse();
	    ucell.G = ucell.GT.Transpose();
	    ucell.lat0 = 1.8897261254578281;
	    for (int i = 0; i < coord.size(); i++)
	    {
	        ucell.atoms[i].label = coord[i].atomname;
	        ucell.atoms[i].na = coord[i].coordinate.size();
	        ucell.atoms[i].tau = new ModuleBase::Vector3<double>[ucell.atoms[i].na];
	        ucell.atoms[i].taud = new ModuleBase::Vector3<double>[ucell.atoms[i].na];
	        for (int j = 0; j < ucell.atoms[i].na; j++)
	        {
	            std::vector<double> this_atom = coord[i].coordinate[j];
	            ucell.atoms[i].tau[j] = ModuleBase::Vector3<double>(this_atom[0], this_atom[1], this_atom[2]);
	            ModuleBase::Mathzone::Cartesian_to_Direct(ucell.atoms[i].tau[j].x,
	                                                      ucell.atoms[i].tau[j].y,
	                                                      ucell.atoms[i].tau[j].z,
	                                                      ucell.a1.x,
	                                                      ucell.a1.y,
	                                                      ucell.a1.z,
	                                                      ucell.a2.x,
	                                                      ucell.a2.y,
	                                                      ucell.a2.z,
	                                                      ucell.a3.x,
	                                                      ucell.a3.y,
	                                                      ucell.a3.z,
	                                                      ucell.atoms[i].taud[j].x,
	                                                      ucell.atoms[i].taud[j].y,
	                                                      ucell.atoms[i].taud[j].z);
	        }
	        ucell.nat += ucell.atoms[i].na;
	    }
	}
	//clear ucell
	void ClearUcell()
	{
	    for (int i = 0; i < ucell.ntype; i++)
	    {
	        delete[] ucell.atoms[i].tau;
	        delete[] ucell.atoms[i].taud;
	    }
	    delete[] ucell.atoms;
	}
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
	//just to set ucell info here, however it is used in the following tests
	GlobalC::ucell.latvec.e11 = 10.0; GlobalC::ucell.latvec.e12 = 0.0; GlobalC::ucell.latvec.e13 = 0.0;
	GlobalC::ucell.latvec.e21 = 0.0; GlobalC::ucell.latvec.e22 = 10.0; GlobalC::ucell.latvec.e23 = 0.0;
	GlobalC::ucell.latvec.e31 = 0.0; GlobalC::ucell.latvec.e32 = 0.0; GlobalC::ucell.latvec.e33 = 10.0;
	GlobalC::ucell.GT = GlobalC::ucell.latvec.Inverse();
	GlobalC::ucell.G = GlobalC::ucell.GT.Transpose();
	GlobalC::ucell.lat0 = 1.8897261254578281;
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
	GlobalV::KSPACING[0] = 0.052918; // 0.52918/Bohr = 1/A
	GlobalV::KSPACING[1] = 0.052918; // 0.52918/Bohr = 1/A
	GlobalV::KSPACING[2] = 0.052918; // 0.52918/Bohr = 1/A
	std::string k_file = "./support/KPT3";
	kv->read_kpoints(k_file);
	EXPECT_EQ(kv->nkstot,343);
	GlobalV::KSPACING[0]=0.0;
	GlobalV::KSPACING[1]=0.0;
	GlobalV::KSPACING[2]=0.0;
}

TEST_F(KlistTest, ReadKpointsKspacing3values)
{
	kv->nspin = 1;
	GlobalV::KSPACING[0] = 0.052918; // 0.52918/Bohr = 1/A
	GlobalV::KSPACING[1] = 0.06; // 0.52918/Bohr = 1/A
	GlobalV::KSPACING[2] = 0.07; // 0.52918/Bohr = 1/A
	std::string k_file = "./support/KPT3";
	kv->read_kpoints(k_file);
	EXPECT_EQ(kv->nkstot,210);
	GlobalV::KSPACING[0]=0.0;
	GlobalV::KSPACING[1]=0.0;
	GlobalV::KSPACING[2]=0.0;
}

TEST_F(KlistTest, ReadKpointsInvalidKspacing3values)
{
	kv->nspin = 1;
	GlobalV::KSPACING[0] = 0.052918; // 0.52918/Bohr = 1/A
	GlobalV::KSPACING[1] = 0; // 0.52918/Bohr = 1/A
	GlobalV::KSPACING[2] = 0.07; // 0.52918/Bohr = 1/A
	std::string k_file = "./support/KPT3";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(kv->read_kpoints(k_file),::testing::ExitedWithCode(0),"");
	output = testing::internal::GetCapturedStdout();
	GlobalV::KSPACING[0]=0.0;
	GlobalV::KSPACING[1]=0.0;
	GlobalV::KSPACING[2]=0.0;
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
	ifs.open("klist_tmp_warning_1");
	std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
	EXPECT_THAT(str, testing::HasSubstr("Can't find File name : arbitrary_1"));
	ifs.close();
	remove("klist_tmp_warning_1");
}

TEST_F(KlistTest, ReadKpointsWarning2)
{
	std::string k_file = "arbitrary_2";
	ofs.open(k_file.c_str());
	ofs<<"ARBITRARY";
	ofs.close();
	kv->nspin = 1;
	GlobalV::ofs_warning.open("klist_tmp_warning_2");
	EXPECT_NO_THROW(kv->read_kpoints(k_file));
	GlobalV::ofs_warning.close();
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
	ofs.open(k_file.c_str());
	ofs<<"KPOINTS"<<std::endl;
	ofs<<"100001"<<std::endl;
	ofs.close();
	kv->nspin = 1;
	GlobalV::ofs_warning.open("klist_tmp_warning_3");
	EXPECT_NO_THROW(kv->read_kpoints(k_file));
	GlobalV::ofs_warning.close();
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
	ofs.open(k_file.c_str());
	ofs<<"KPOINTS"<<std::endl;
	ofs<<"0"<<std::endl;
	ofs<<"arbitrary"<<std::endl;
	ofs.close();
	kv->nspin = 1;
	GlobalV::ofs_warning.open("klist_tmp_warning_4");
	EXPECT_NO_THROW(kv->read_kpoints(k_file));
	GlobalV::ofs_warning.close();
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
	ofs.open(k_file.c_str());
	ofs<<"KPOINTS"<<std::endl;
	ofs<<"100000"<<std::endl;
	ofs<<"Line_Direct"<<std::endl;
	ofs.close();
	kv->nspin = 1;
	ModuleSymmetry::Symmetry::symm_flag = 1;
	GlobalV::ofs_warning.open("klist_tmp_warning_7");
	EXPECT_NO_THROW(kv->read_kpoints(k_file));
	GlobalV::ofs_warning.close();
	ifs.open("klist_tmp_warning_7");
	std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
	EXPECT_THAT(str, testing::HasSubstr("Line mode of k-points is open, please set symmetry to 0 or -1"));
	ifs.close();
	remove("klist_tmp_warning_7");
	remove("arbitrary_7");
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
	kv->nkstot = 1;
	kv->nks = 2;
	kv->renew(kv->nkstot);
	kv->kvec_c[0].x = 0;
	kv->kvec_c[0].y = 0;
	kv->kvec_c[0].z = 0;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(kv->print_klists(GlobalV::ofs_running),
			::testing::ExitedWithCode(0),"");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("nkstot < nks"));
}

TEST_F(KlistTest, SetBothKvecFinalSCF)
{
	kv->nspin = 1;
	kv->nkstot = 1;
	kv->nks = 1;
	kv->renew(kv->nkstot);
	kv->kvec_d[0].x = 0.0;
	kv->kvec_d[0].y = 0.0;
	kv->kvec_d[0].z = 0.0;
	kv->kvec_c[0].x = 0.0;
	kv->kvec_c[0].y = 0.0;
	kv->kvec_c[0].z = 0.0;
	std::string skpt;
	GlobalV::FINAL_SCF = 1;
	kv->kd_done = 0;
	kv->kc_done = 0;
	// case 1
	kv->k_nkstot = 0;
	kv->set_both_kvec(GlobalC::ucell.G,GlobalC::ucell.latvec,skpt);
	EXPECT_TRUE(kv->kd_done);
	EXPECT_TRUE(kv->kc_done);
	// case 2
	kv->k_nkstot = 1;
	kv->k_kword = "D";
	kv->set_both_kvec(GlobalC::ucell.G,GlobalC::ucell.latvec,skpt);
	EXPECT_TRUE(kv->kd_done);
	EXPECT_TRUE(kv->kc_done);
	// case 3
	kv->k_kword = "C";
	kv->set_both_kvec(GlobalC::ucell.G,GlobalC::ucell.latvec,skpt);
	EXPECT_TRUE(kv->kc_done);
	EXPECT_TRUE(kv->kd_done);
	// case 4
	GlobalV::ofs_warning.open("klist_tmp_warning_8");
	kv->k_kword = "arbitrary";
	kv->set_both_kvec(GlobalC::ucell.G,GlobalC::ucell.latvec,skpt);
	GlobalV::ofs_warning.close();
	ifs.open("klist_tmp_warning_8");
	std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
	EXPECT_THAT(str, testing::HasSubstr("Error : neither Cartesian nor Direct kpoint."));
	ifs.close();
	remove("klist_tmp_warning_8");
}

TEST_F(KlistTest, SetBothKvec)
{
	kv->nspin = 1;
	kv->nkstot = 1;
	kv->nks = 1;
	kv->renew(kv->nkstot);
	kv->kvec_d[0].x = 0.0;
	kv->kvec_d[0].y = 0.0;
	kv->kvec_d[0].z = 0.0;
	kv->kc_done = 0;
	kv->kd_done = 1;
	std::string skpt;
	GlobalV::FINAL_SCF = 0;
	kv->set_both_kvec(GlobalC::ucell.G,GlobalC::ucell.latvec,skpt);
	EXPECT_TRUE(kv->kc_done);
	kv->kc_done = 1;
	kv->kd_done = 0;
	kv->set_both_kvec(GlobalC::ucell.G,GlobalC::ucell.latvec,skpt);
	EXPECT_TRUE(kv->kd_done);
}

TEST_F(KlistTest, NormalizeWk)
{
	kv->nspin = 1;
	kv->nkstot = 2;
	kv->nks = 2;
	kv->renew(kv->nkstot);
	kv->wk[0] = 1.0;
	kv->wk[1] = 1.0;
	int deg = 2;
	kv->normalize_wk(deg);
	EXPECT_DOUBLE_EQ(kv->wk[0],1.0);
	EXPECT_DOUBLE_EQ(kv->wk[1],1.0);
}

TEST_F(KlistTest, UpdateUseIBZ)
{
	kv->nspin = 1;
	kv->nkstot = 3;
	kv->nks = 3;
	kv->renew(kv->nkstot);
	kv->nkstot_ibz = 2;
	kv->kvec_d_ibz.resize(kv->nkstot_ibz);
	kv->wk_ibz.resize(kv->nkstot_ibz);
	kv->update_use_ibz();
	EXPECT_EQ(kv->nkstot,2);
	EXPECT_EQ(kv->kvec_d.size(),2);
	EXPECT_TRUE(kv->kd_done);
	EXPECT_FALSE(kv->kc_done);
}

TEST_F(KlistTest, IbzKpoint)
{
	//construct cell and symmetry
	ModuleSymmetry::Symmetry symm;
	construct_ucell(stru_lib[0]);
	GlobalV::ofs_running.open("tmp_klist_3");
	symm.analy_sys(ucell,GlobalV::ofs_running);
	//read KPT
	std::string k_file = "./support/KPT1";
	kv->nspin = 1;
	kv->read_kpoints(k_file);
	EXPECT_EQ(kv->nkstot,512);
	//calculate ibz_kpoint
	std::string skpt;
    ModuleSymmetry::Symmetry::symm_flag = 1;
    bool match = true;
    kv->ibz_kpoint(symm, ModuleSymmetry::Symmetry::symm_flag, skpt, ucell, match);
	EXPECT_EQ(kv->nkstot_ibz,35);
	GlobalV::ofs_running<<skpt<<std::endl;
	GlobalV::ofs_running.close();
	ClearUcell();
	remove("tmp_klist_3");
}

TEST_F(KlistTest, IbzKpointIsMP)
{
	//construct cell and symmetry
	ModuleSymmetry::Symmetry symm;
	construct_ucell(stru_lib[0]);
	GlobalV::ofs_running.open("tmp_klist_4");
	symm.analy_sys(ucell,GlobalV::ofs_running);
	//read KPT
	std::string k_file = "./support/KPT1";
	kv->nspin = 1;
	kv->read_kpoints(k_file);
	EXPECT_EQ(kv->nkstot,512);
	EXPECT_TRUE(kv->is_mp);
	//calculate ibz_kpoint
	std::string skpt;
    ModuleSymmetry::Symmetry::symm_flag = 0;
    bool match = true;
    kv->ibz_kpoint(symm, ModuleSymmetry::Symmetry::symm_flag, skpt, ucell, match);
	EXPECT_EQ(kv->nkstot_ibz,260);
	GlobalV::ofs_running<<skpt<<std::endl;
	GlobalV::ofs_running.close();
	ClearUcell();
	remove("tmp_klist_4");
}

#undef private

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	testing::InitGoogleTest(&argc, argv);

	int nproc;
	int my_rank;
	
	MPI_Comm_size(MPI_COMM_WORLD,&nproc);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	
	int result = RUN_ALL_TESTS();
	
	MPI_Finalize();
	
	return result;
}