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
#include "module_cell/pseudo.h"
#include "module_cell/setup_nonlocal.h"
#include "module_hamilt_pw/hamilt_pwdft/parallel_grid.h"
#include "module_cell/parallel_kpoints.h"
#include "module_io/berryphase.h"
#include "module_basis/module_ao/ORB_gaunt_table.h"

bool berryphase::berry_phase_flag=0;

pseudo::pseudo(){}
pseudo::~pseudo(){}
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
Soc::~Soc()
{
}
Fcoef::~Fcoef()
{
}

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
 *   - Set
 *     - this is a "kind of" integerated test
 *       for set() and mpi_k()
 *   - SetAfterVC
 *     - this is a "kind of" integerated test
 *       for set_after_vc() and mpi_k_after_vc()
 *     - a bug is found from here, that is,
 *       KPAR > 1 is not support yet in vc-relax calculation
 *       due to the size of kvec_d, kvec_c being nks, rather
 *       than nkstot in set_both_kvec_after_vc
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


class KlistParaTest : public testing::Test
{
protected:
	std::unique_ptr<K_Vectors> kv{new K_Vectors};
	std::ifstream ifs;
	std::ofstream ofs;
	std::ofstream ofs_running;
	std::string output;

	//used to construct cell and analyse its symmetry
	void construct_ucell(stru_ &stru)
	{
	    std::vector<atomtype_> coord = stru.all_type;
	    GlobalC::ucell.a1 = ModuleBase::Vector3<double>(stru.cell[0], stru.cell[1], stru.cell[2]);
	    GlobalC::ucell.a2 = ModuleBase::Vector3<double>(stru.cell[3], stru.cell[4], stru.cell[5]);
	    GlobalC::ucell.a3 = ModuleBase::Vector3<double>(stru.cell[6], stru.cell[7], stru.cell[8]);
	    GlobalC::ucell.ntype = stru.all_type.size();
	    GlobalC::ucell.atoms = new Atom[GlobalC::ucell.ntype];
	    GlobalC::ucell.nat = 0;
	    GlobalC::ucell.latvec.e11 = GlobalC::ucell.a1.x; GlobalC::ucell.latvec.e12 = GlobalC::ucell.a1.y; GlobalC::ucell.latvec.e13 = GlobalC::ucell.a1.z;
	    GlobalC::ucell.latvec.e21 = GlobalC::ucell.a2.x; GlobalC::ucell.latvec.e22 = GlobalC::ucell.a2.y; GlobalC::ucell.latvec.e23 = GlobalC::ucell.a2.z;
	    GlobalC::ucell.latvec.e31 = GlobalC::ucell.a3.x; GlobalC::ucell.latvec.e32 = GlobalC::ucell.a3.y; GlobalC::ucell.latvec.e33 = GlobalC::ucell.a3.z;
	    GlobalC::ucell.GT = GlobalC::ucell.latvec.Inverse();
	    GlobalC::ucell.G = GlobalC::ucell.GT.Transpose();
	    GlobalC::ucell.lat0 = 1.8897261254578281;
	    for (int i = 0; i < coord.size(); i++)
	    {
	        GlobalC::ucell.atoms[i].label = coord[i].atomname;
	        GlobalC::ucell.atoms[i].na = coord[i].coordinate.size();
	        GlobalC::ucell.atoms[i].tau = new ModuleBase::Vector3<double>[GlobalC::ucell.atoms[i].na];
	        GlobalC::ucell.atoms[i].taud = new ModuleBase::Vector3<double>[GlobalC::ucell.atoms[i].na];
	        for (int j = 0; j < GlobalC::ucell.atoms[i].na; j++)
	        {
	            std::vector<double> this_atom = coord[i].coordinate[j];
	            GlobalC::ucell.atoms[i].tau[j] = ModuleBase::Vector3<double>(this_atom[0], this_atom[1], this_atom[2]);
	            ModuleBase::Mathzone::Cartesian_to_Direct(GlobalC::ucell.atoms[i].tau[j].x,
	                                                      GlobalC::ucell.atoms[i].tau[j].y,
	                                                      GlobalC::ucell.atoms[i].tau[j].z,
	                                                      GlobalC::ucell.a1.x,
	                                                      GlobalC::ucell.a1.y,
	                                                      GlobalC::ucell.a1.z,
	                                                      GlobalC::ucell.a2.x,
	                                                      GlobalC::ucell.a2.y,
	                                                      GlobalC::ucell.a2.z,
	                                                      GlobalC::ucell.a3.x,
	                                                      GlobalC::ucell.a3.y,
	                                                      GlobalC::ucell.a3.z,
	                                                      GlobalC::ucell.atoms[i].taud[j].x,
	                                                      GlobalC::ucell.atoms[i].taud[j].y,
	                                                      GlobalC::ucell.atoms[i].taud[j].z);
	        }
	        GlobalC::ucell.nat += GlobalC::ucell.atoms[i].na;
	    }
	}
	//clear GlobalC::ucell
	void ClearUcell()
	{
	    for (int i = 0; i < GlobalC::ucell.ntype; i++)
	    {
	        delete[] GlobalC::ucell.atoms[i].tau;
	        delete[] GlobalC::ucell.atoms[i].taud;
	    }
	    delete[] GlobalC::ucell.atoms;
	}
};

#ifdef __MPI
TEST_F(KlistParaTest,Set)
{
	//construct cell and symmetry
	ModuleSymmetry::Symmetry symm;
	construct_ucell(stru_lib[0]);
    if (GlobalV::MY_RANK == 0) GlobalV::ofs_running.open("tmp_klist_5");
    symm.analy_sys(GlobalC::ucell.lat, GlobalC::ucell.st, GlobalC::ucell.atoms, GlobalV::ofs_running);
	//read KPT
	std::string k_file = "./support/KPT1";
	//set klist
	kv->nspin = 1;
	GlobalV::NSPIN = 1;
	if(GlobalV::NPROC == 4) {GlobalV::KPAR = 2;}
	Parallel_Global::init_pools();
	ModuleSymmetry::Symmetry::symm_flag=1;
	kv->set(symm,k_file,kv->nspin,GlobalC::ucell.G,GlobalC::ucell.latvec);
	EXPECT_EQ(kv->nkstot,35);
	EXPECT_TRUE(kv->kc_done);
	EXPECT_TRUE(kv->kd_done);
	if(GlobalV::NPROC == 4)
	{
		if(GlobalV::MY_RANK==0) EXPECT_EQ(kv->nks,18);
		if(GlobalV::MY_RANK==1) EXPECT_EQ(kv->nks,18);
		if(GlobalV::MY_RANK==2) EXPECT_EQ(kv->nks,17);
		if(GlobalV::MY_RANK==3) EXPECT_EQ(kv->nks,17);
	}
	ClearUcell();
	if(GlobalV::MY_RANK==0)
	{
		GlobalV::ofs_running.close();
		remove("tmp_klist_5");
		remove("kpoints");
	}
}

TEST_F(KlistParaTest,SetAfterVC)
{
	//construct cell and symmetry
	ModuleSymmetry::Symmetry symm;
	construct_ucell(stru_lib[0]);
	if(GlobalV::MY_RANK == 0) GlobalV::ofs_running.open("tmp_klist_6");
    symm.analy_sys(GlobalC::ucell.lat, GlobalC::ucell.st, GlobalC::ucell.atoms, GlobalV::ofs_running);
	//read KPT
	std::string k_file = "./support/KPT1";
	//set klist
	kv->nspin = 1;
	GlobalV::NSPIN = 1;
	if(GlobalV::NPROC == 4) {GlobalV::KPAR = 1;}
	Parallel_Global::init_pools();
	ModuleSymmetry::Symmetry::symm_flag=1;
	kv->set(symm,k_file,kv->nspin,GlobalC::ucell.G,GlobalC::ucell.latvec);
	EXPECT_EQ(kv->nkstot,35);
	EXPECT_TRUE(kv->kc_done);
	EXPECT_TRUE(kv->kd_done);
	if(GlobalV::NPROC == 4)
	{
		if(GlobalV::MY_RANK==0) EXPECT_EQ(kv->nks,35);
		if(GlobalV::MY_RANK==1) EXPECT_EQ(kv->nks,35);
		if(GlobalV::MY_RANK==2) EXPECT_EQ(kv->nks,35);
		if(GlobalV::MY_RANK==3) EXPECT_EQ(kv->nks,35);
	}
	//call set_after_vc here
	kv->kc_done = 0;
	kv->set_after_vc(symm,k_file,kv->nspin,GlobalC::ucell.G,GlobalC::ucell.latvec);
	EXPECT_TRUE(kv->kc_done);
	EXPECT_TRUE(kv->kd_done);
	//clear
	ClearUcell();
	if(GlobalV::MY_RANK==0)
	{
		GlobalV::ofs_running.close();
		remove("tmp_klist_6");
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
