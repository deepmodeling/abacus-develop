#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "memory"
#include "module_base/mathzone.h"
#include "module_base/global_variable.h"
#include "module_cell/unitcell.h"
#include<vector>
#include<valarray>
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
 *   - Index
 *     - set_iat2iait(): set index relations in two arrays of Unitcell: iat2it[nat], iat2ia[nat]
 *     - iat2iait(): depends on the above function, but can find both ia & it from iat
 *     - ijat2iaitjajt(): find ia, it, ja, jt from ijat (ijat_max = nat*nat)
 *         which collapses it, ia, jt, ja loop into a single loop
 *     - step_ia(): periodically set ia to 0 when ia reaches atom[it].na - 1
 *     - step_it(): periodically set it to 0 when it reaches ntype -1 
 *     - step_iait(): return true only the above two conditions are true
 *     - step_jajtiait(): return ture only two of the above function (for i and j) are true
 *   - CheckDTau
 *     - check_dtau(): move all atomic coordinates into the first unitcell, i.e. in between [0,1)
 *   - CheckTau
 *     - check_tau(): check if any "two atoms are too close"
 *   - SelectiveDynamics
 *     - if_atoms_can_move():it is true if any coordinates of any atom can move, i.e. mbl = 1
 *   - SaveCartesianPosition
 *     - save_cartesian_position(): make a copy of atomic Cartesian coordinates: tau
 *   - SaveCartesianPositionOriginal
 *     - save_cartesian_position_original(): make a copy of original atomic Cartesian coordinates: tau_original
 *       tau_original means without periodic adjustment
 *   - PeriodicBoundaryAdjustment
 *     - periodic_boundary_adjustment(): move atoms inside the unitcell after relaxation
 *   - PrintCell
 *     - print_cell(ofs): print basic cell info into ofs
 *   - PrintCellCif
 *     - print_cell_cif(fn): print cif file (named as fn)
 *   - PrintSTRU
 *     - print_stru_file(): print STRU file of ABACUS
 *   - PrintTau
 *     - print_tau(): print atomic coordinates, magmom and initial velocities
 *   - PrintUnitcellPseudo
 *     - Actually an integrated function to call UnitCell::print_cell and Atom::print_Atom
 *   - UpdatePosTau1
 *     - update_pos_tau(const double* pos)
 *   - UpdatePosTau2
 *     - update_pos_tau(const ModuleBase::Vector3<double>* posd_in)
 *   - UpdateVel
 *     - update_vel(const ModuleBase::Vector3<double>* vel_in)
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

class UcellTest : public ::testing::TestWithParam<UcellTestPrepare>{};

TEST_P(UcellTest,Index)
{
	UcellTestPrepare utp = GetParam();
	std::unique_ptr<UnitCell> ucell;
	GlobalV::relax_new = utp.relax_new;
	ucell = utp.SetUcellInfo();
	if(utp.system_name == "C1H2-Index")
	{
		//test set_iat2itia
		ucell->set_iat2itia();
		int iat = 0;
		for (int it=0; it<utp.natom.size(); ++it)
		{
			for(int ia=0; ia<utp.natom[it];++ia)
			{
				EXPECT_EQ(ucell->iat2it[iat],it);
				EXPECT_EQ(ucell->iat2ia[iat],ia);
				// test iat2iait
				int ia_beg, it_beg;
				ucell->iat2iait(iat,&ia_beg,&it_beg);
				EXPECT_EQ(it_beg,it);
				EXPECT_EQ(ia_beg,ia);
				++iat;
			}
		}
		// test iat2iait: case of (iat >= nat)
		int ia_beg2;
		int it_beg2;
		long long iat2 = ucell->nat + 1;
		EXPECT_FALSE(ucell->iat2iait(iat2,&ia_beg2,&it_beg2));
		// test ijat2iaitjajt, step_jajtiait, step_iat, step_ia, step_it
		int ia_test;
		int it_test;
		int ja_test;
		int jt_test;
		int ia_test2 = 0;
		int it_test2 = 0;
		int ja_test2 = 0;
		int jt_test2 = 0;
		long long ijat = 0;
		for(int it=0; it<utp.natom.size(); ++it)
		{
			for(int ia=0; ia<utp.natom[it];++ia)
			{
				for(int jt=0; jt<utp.natom.size(); ++jt)
				{
					for(int ja=0; ja<utp.natom[jt];++ja)
					{
						ucell->ijat2iaitjajt(ijat,&ia_test,&it_test,&ja_test,&jt_test);
						EXPECT_EQ(ia_test,ia);
						EXPECT_EQ(it_test,it);
						EXPECT_EQ(ja_test,ja);
						EXPECT_EQ(jt_test,jt);
						++ijat;
						if (it_test == utp.natom.size()-1 &&
							ia_test == utp.natom[it]-1 &&
							jt_test == utp.natom.size()-1 &&
							ja_test == utp.natom[jt]-1)
						{
							EXPECT_TRUE(ucell->step_jajtiait(&ja_test,&jt_test,&ia_test,&it_test));
						}
						else
						{
							EXPECT_FALSE(ucell->step_jajtiait(&ja_test,&jt_test,&ia_test,&it_test));
						}
					}
				}
			}
		}
	}
}

TEST_P(UcellTest,CheckDTau)
{
	UcellTestPrepare utp = GetParam();
	std::unique_ptr<UnitCell> ucell;
	GlobalV::relax_new = utp.relax_new;
	ucell = utp.SetUcellInfo();
	if(utp.system_name == "C1H2-CheckDTau")
	{
		ucell->check_dtau();
		for (int it=0; it<utp.natom.size(); ++it)
		{
			for(int ia=0; ia<utp.natom[it];++ia)
			{
				EXPECT_GE(ucell->atoms[it].taud[ia].x,0);
				EXPECT_GE(ucell->atoms[it].taud[ia].y,0);
				EXPECT_GE(ucell->atoms[it].taud[ia].z,0);
				EXPECT_LT(ucell->atoms[it].taud[ia].x,1);
				EXPECT_LT(ucell->atoms[it].taud[ia].y,1);
				EXPECT_LT(ucell->atoms[it].taud[ia].z,1);
			}
		}
	}
}

TEST_P(UcellTest,CheckTau)
{
	UcellTestPrepare utp = GetParam();
	std::unique_ptr<UnitCell> ucell;
	GlobalV::relax_new = utp.relax_new;
	ucell = utp.SetUcellInfo();
	if(utp.system_name == "C1H2-CheckTau")
	{
		GlobalV::ofs_warning.open("checktau_warning");
		ucell->check_tau();
		GlobalV::ofs_warning.close();
		std::ifstream ifs;
		ifs.open("checktau_warning");
    		std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
    		EXPECT_THAT(str, testing::HasSubstr("two atoms are too close!"));
		ifs.close();
		remove("checktau_warning");
	}
}

TEST_P(UcellTest,SelectiveDynamics)
{
	UcellTestPrepare utp = GetParam();
	std::unique_ptr<UnitCell> ucell;
	GlobalV::relax_new = utp.relax_new;
	ucell = utp.SetUcellInfo();
	if(utp.system_name == "C1H2-SD")
	{
		EXPECT_TRUE(ucell->if_atoms_can_move());
	}
}

TEST_P(UcellTest,SaveCartesianPosition)
{
	UcellTestPrepare utp = GetParam();
	std::unique_ptr<UnitCell> ucell;
	GlobalV::relax_new = utp.relax_new;
	ucell = utp.SetUcellInfo();
	//the first realization
	double* pos = new double[3*ucell->nat];
	ucell->save_cartesian_position(pos);
	//another realization
	ModuleBase::Vector3<double>* pos1 = new ModuleBase::Vector3<double>[ucell->nat];
	ucell->save_cartesian_position(pos1);
	int iat = 0;
	for(int it=0; it<utp.natom.size(); ++it)
	{
		for(int ia=0; ia<utp.natom[it];++ia)
		{
			//first
			EXPECT_DOUBLE_EQ(pos[3*iat],ucell->atoms[it].tau[ia].x*ucell->lat0);
			EXPECT_DOUBLE_EQ(pos[3*iat+1],ucell->atoms[it].tau[ia].y*ucell->lat0);
			EXPECT_DOUBLE_EQ(pos[3*iat+2],ucell->atoms[it].tau[ia].z*ucell->lat0);
			//second
			EXPECT_DOUBLE_EQ(pos1[iat].x,ucell->atoms[it].tau[ia].x*ucell->lat0);
			EXPECT_DOUBLE_EQ(pos1[iat].y,ucell->atoms[it].tau[ia].y*ucell->lat0);
			EXPECT_DOUBLE_EQ(pos1[iat].z,ucell->atoms[it].tau[ia].z*ucell->lat0);
			++iat;
		}
	}
	delete[] pos;
	delete[] pos1;
}

TEST_P(UcellTest,SaveCartesianPositionOriginal)
{
	UcellTestPrepare utp = GetParam();
	std::unique_ptr<UnitCell> ucell;
	GlobalV::relax_new = utp.relax_new;
	ucell = utp.SetUcellInfo();
	//the first realization
	double* pos = new double[3*ucell->nat];
	ucell->save_cartesian_position_original(pos);
	//another realization
	ModuleBase::Vector3<double>* pos1 = new ModuleBase::Vector3<double>[ucell->nat];
	ucell->save_cartesian_position_original(pos1);
	int iat = 0;
	for(int it=0; it<utp.natom.size(); ++it)
	{
		for(int ia=0; ia<utp.natom[it];++ia)
		{
			//first
			EXPECT_DOUBLE_EQ(pos[3*iat],ucell->atoms[it].tau_original[ia].x*ucell->lat0);
			EXPECT_DOUBLE_EQ(pos[3*iat+1],ucell->atoms[it].tau_original[ia].y*ucell->lat0);
			EXPECT_DOUBLE_EQ(pos[3*iat+2],ucell->atoms[it].tau[ia].z*ucell->lat0);
			//second
			EXPECT_DOUBLE_EQ(pos1[iat].x,ucell->atoms[it].tau_original[ia].x*ucell->lat0);
			EXPECT_DOUBLE_EQ(pos1[iat].y,ucell->atoms[it].tau_original[ia].y*ucell->lat0);
			EXPECT_DOUBLE_EQ(pos1[iat].z,ucell->atoms[it].tau_original[ia].z*ucell->lat0);
			++iat;
		}
	}
	delete[] pos;
	delete[] pos1;
}

TEST_P(UcellTest,PeriodicBoundaryAdjustment)
{
	UcellTestPrepare utp = GetParam();
	std::unique_ptr<UnitCell> ucell;
	GlobalV::relax_new = utp.relax_new;
	ucell = utp.SetUcellInfo();
	if(utp.system_name == "C1H2-PBA")
	{
		testing::internal::CaptureStdout();
		EXPECT_EXIT(ucell->periodic_boundary_adjustment(),
				::testing::ExitedWithCode(0),"");
		std::string output = testing::internal::GetCapturedStdout();
		EXPECT_THAT(output,testing::HasSubstr("the movement of atom is larger than the length of cell"));
	}
	else if (utp.system_name == "C1H2-Index")
	{
		EXPECT_NO_THROW(ucell->periodic_boundary_adjustment());
	}
}

TEST_P(UcellTest,PrintCell)
{
	UcellTestPrepare utp = GetParam();
	std::unique_ptr<UnitCell> ucell;
	GlobalV::relax_new = utp.relax_new;
	ucell = utp.SetUcellInfo();
	if(utp.system_name == "C1H2-Index")
	{
		std::ofstream ofs;
		ofs.open("printcell.log");
		ucell->print_cell(ofs);
		ofs.close();
		std::ifstream ifs;
    		ifs.open("printcell.log");
    		std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
    		EXPECT_THAT(str, testing::HasSubstr("latName = bcc"));
    		EXPECT_THAT(str, testing::HasSubstr("ntype = 2"));
    		EXPECT_THAT(str, testing::HasSubstr("nat = 3"));
    		EXPECT_THAT(str, testing::HasSubstr("GGT :"));
    		EXPECT_THAT(str, testing::HasSubstr("omega = 6748.33"));
		remove("printcell.log");
	}
}

TEST_P(UcellTest,PrintUnitcellPseudo)
{
	UcellTestPrepare utp = GetParam();
	std::unique_ptr<UnitCell> ucell;
	GlobalV::relax_new = utp.relax_new;
	ucell = utp.SetUcellInfo();
	if(utp.system_name == "C1H2-Index")
	{
		GlobalV::test_pseudo_cell = 1;
		std::string fn="printcell.log";
		ucell->print_unitcell_pseudo(fn);
		std::ifstream ifs;
    		ifs.open("printcell.log");
    		std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
    		EXPECT_THAT(str, testing::HasSubstr("latName = bcc"));
    		EXPECT_THAT(str, testing::HasSubstr("ntype = 2"));
    		EXPECT_THAT(str, testing::HasSubstr("nat = 3"));
    		EXPECT_THAT(str, testing::HasSubstr("GGT :"));
    		EXPECT_THAT(str, testing::HasSubstr("omega = 6748.33"));
    		EXPECT_THAT(str, testing::HasSubstr("label = C"));
    		EXPECT_THAT(str, testing::HasSubstr("mass = 12"));
    		EXPECT_THAT(str, testing::HasSubstr("atom_position(cartesian) Dimension = 1"));
    		EXPECT_THAT(str, testing::HasSubstr("label = H"));
    		EXPECT_THAT(str, testing::HasSubstr("mass = 1"));
    		EXPECT_THAT(str, testing::HasSubstr("atom_position(cartesian) Dimension = 2"));
		remove("printcell.log");
	}
}

TEST_P(UcellTest,PrintCellCif)
{
	UcellTestPrepare utp = GetParam();
	std::unique_ptr<UnitCell> ucell;
	GlobalV::relax_new = utp.relax_new;
	ucell = utp.SetUcellInfo();
	if(utp.system_name == "C1H2-Index")
	{
		GlobalV::test_unitcell = 1;
		GlobalV::global_out_dir = "./";
		std::string fn = "printcell.cif";
		ucell->print_cell_cif(fn);
		std::ifstream ifs;
    		ifs.open("printcell.cif");
    		std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
    		EXPECT_THAT(str, testing::HasSubstr("data_bcc"));
    		EXPECT_THAT(str, testing::HasSubstr("_audit_creation_method generated by ABACUS"));
		remove("printcell.cif");
	}
}

TEST_P(UcellTest,PrintSTRU)
{
	UcellTestPrepare utp = GetParam();
	std::unique_ptr<UnitCell> ucell;
	GlobalV::relax_new = utp.relax_new;
	ucell = utp.SetUcellInfo();
	if(utp.system_name == "C1H2-Index")
	{
		//Cartesian type of coordinates
		std::string fn = "C1H2_STRU";
		int type = 1; // for Cartesian
		int level = 1; //print velocity in STRU
		ucell->print_stru_file(fn,type,level);
		std::ifstream ifs;
		ifs.open("C1H2_STRU");
    		std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
    		EXPECT_THAT(str, testing::HasSubstr("C 12 C.upf upf201"));
    		EXPECT_THAT(str, testing::HasSubstr("H 1 H.upf upf201"));
    		EXPECT_THAT(str, testing::HasSubstr("Cartesian"));
    		EXPECT_THAT(str, testing::HasSubstr("H #label"));
    		EXPECT_THAT(str, testing::HasSubstr("0 #magnetism"));
    		EXPECT_THAT(str, testing::HasSubstr("2 #number of atoms"));
    		EXPECT_THAT(str, testing::HasSubstr("1.5  1.5  1.5  m  0  0  0  v  0.1  0.1  0.1"));
    		EXPECT_THAT(str, testing::HasSubstr("0.5  0.5  0.5  m  0  0  1  v  0.1  0.1  0.1"));
		str.clear();
		ifs.close();
		remove("C1H2_STRU");
		//direct type of coordinates
		type = 2; //for direct
		ucell->print_stru_file(fn,type,level);
		ifs.open("C1H2_STRU");
		str= {(std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>()};
    		EXPECT_THAT(str, testing::HasSubstr("C 12 C.upf upf201"));
    		EXPECT_THAT(str, testing::HasSubstr("H 1 H.upf upf201"));
    		EXPECT_THAT(str, testing::HasSubstr("Direct"));
    		EXPECT_THAT(str, testing::HasSubstr("H #label"));
    		EXPECT_THAT(str, testing::HasSubstr("0 #magnetism"));
    		EXPECT_THAT(str, testing::HasSubstr("2 #number of atoms"));
    		EXPECT_THAT(str, testing::HasSubstr("0.15  0.15  0.15  m  0  0  0  v  0.1  0.1  0.1"));
    		EXPECT_THAT(str, testing::HasSubstr("0.05  0.05  0.05  m  0  0  1  v  0.1  0.1  0.1"));
		ifs.close();
		remove("C1H2_STRU");
	}
}

TEST_P(UcellTest,PrintTau)
{
	UcellTestPrepare utp = GetParam();
	std::unique_ptr<UnitCell> ucell;
	GlobalV::relax_new = utp.relax_new;
	ucell = utp.SetUcellInfo();
	if(utp.system_name == "C1H2-Index")
	{
		GlobalV::ofs_running.open("print_tau_direct");
		EXPECT_EQ(ucell->Coordinate,"Direct");
		ucell->print_tau();
		GlobalV::ofs_running.close();
		std::ifstream ifs;
		ifs.open("print_tau_direct");
    		std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
    		EXPECT_THAT(str, testing::HasSubstr("DIRECT COORDINATES"));
    		EXPECT_THAT(str, testing::HasSubstr("taud_C1                 0.1                 0.1                 0.1"));
		ifs.close();
		remove("print_tau_direct");
	}
	else if(utp.system_name == "C1H2-Cartesian")
	{
		GlobalV::ofs_running.open("print_tau_Cartesian");
		EXPECT_EQ(ucell->Coordinate,"Cartesian");
		ucell->print_tau();
		GlobalV::ofs_running.close();
		std::ifstream ifs;
		ifs.open("print_tau_Cartesian");
    		std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
    		EXPECT_THAT(str, testing::HasSubstr("CARTESIAN COORDINATES"));
    		EXPECT_THAT(str, testing::HasSubstr("tauc_C1                   1                   1                   1"));
		ifs.close();
		remove("print_tau_Cartesian");
	}
}

TEST_P(UcellTest,UpdateVel)
{
	UcellTestPrepare utp = GetParam();
	std::unique_ptr<UnitCell> ucell;
	GlobalV::relax_new = utp.relax_new;
	ucell = utp.SetUcellInfo();
	if(utp.system_name == "C1H2-Index")
	{
		ModuleBase::Vector3<double>* vel_in = new ModuleBase::Vector3<double>[ucell->nat];
		for(int iat=0; iat<ucell->nat; ++iat)
		{
			vel_in[iat].set(iat*0.1,iat*0.1,iat*0.1);
		}
		ucell->update_vel(vel_in);
		for(int iat=0; iat<ucell->nat; ++iat)
		{
			EXPECT_DOUBLE_EQ(vel_in[iat].x,0.1*iat);
			EXPECT_DOUBLE_EQ(vel_in[iat].y,0.1*iat);
			EXPECT_DOUBLE_EQ(vel_in[iat].z,0.1*iat);
		}
		delete[] vel_in;
	}
}

TEST_P(UcellTest,UpdatePosTau1)
{
	UcellTestPrepare utp = GetParam();
	std::unique_ptr<UnitCell> ucell;
	GlobalV::relax_new = utp.relax_new;
	ucell = utp.SetUcellInfo();
	if(utp.system_name == "C1H2-Index")
	{
		double* pos_in = new double[ucell->nat*3];
		ucell->set_iat2itia();
		for(int iat=0; iat<ucell->nat; ++iat)
		{
			pos_in[iat*3] = 0.01*ucell->lat0;
			pos_in[iat*3+1] = 0.01*ucell->lat0;
			pos_in[iat*3+2] = 0.01*ucell->lat0;
		}
		ucell->update_pos_tau(pos_in);
		for(int iat=0; iat<ucell->nat; ++iat)
		{
			int it, ia;
			ucell->iat2iait(iat,&ia,&it);
			if(ucell->atoms[it].mbl[ia].x != 0) EXPECT_DOUBLE_EQ(ucell->atoms[it].tau[ia].x,0.01);
			if(ucell->atoms[it].mbl[ia].y != 0) EXPECT_DOUBLE_EQ(ucell->atoms[it].tau[ia].y,0.01);
			if(ucell->atoms[it].mbl[ia].z != 0) EXPECT_DOUBLE_EQ(ucell->atoms[it].tau[ia].z,0.01);
		}
		delete[] pos_in;
	}
}

TEST_P(UcellTest,UpdatePosTau2)
{
	UcellTestPrepare utp = GetParam();
	std::unique_ptr<UnitCell> ucell;
	GlobalV::relax_new = utp.relax_new;
	ucell = utp.SetUcellInfo();
	if(utp.system_name == "C1H2-Index")
	{
		ModuleBase::Vector3<double>* pos_in = new ModuleBase::Vector3<double>[ucell->nat];
		ucell->set_iat2itia();
		for(int iat=0; iat<ucell->nat; ++iat)
		{
			pos_in[iat].set(0.01*ucell->lat0,0.01*ucell->lat0,0.01*ucell->lat0);
		}
		ucell->update_pos_tau(pos_in);
		for(int iat=0; iat<ucell->nat; ++iat)
		{
			int it, ia;
			ucell->iat2iait(iat,&ia,&it);
			if(ucell->atoms[it].mbl[ia].x != 0) EXPECT_DOUBLE_EQ(ucell->atoms[it].tau[ia].x,0.01);
			if(ucell->atoms[it].mbl[ia].y != 0) EXPECT_DOUBLE_EQ(ucell->atoms[it].tau[ia].y,0.01);
			if(ucell->atoms[it].mbl[ia].z != 0) EXPECT_DOUBLE_EQ(ucell->atoms[it].tau[ia].z,0.01);
		}
		delete[] pos_in;
	}
}

INSTANTIATE_TEST_SUITE_P(UcellForTest,UcellTest,::testing::Values(
			UcellTestPrepare("C1H2-Index",	//system-name
				"bcc",		//latname
				2,		//lmaxmax
				true,		//init_vel
				true,		//selective_dyanmics
				true,		//relax_new
				"volume",	//fixed_axes
				1.8897261254578281, //lat0
				{10.0,0.0,0.0,	//latvec
				 0.0,10.0,0.0,
				 0.0,0.0,10.0},
				{"C","H"},	//elements
				{"C.upf","H.upf"},	//upf file
				{"upf201","upf201"},	//upf types
				{"C.orb","H.orb"},	//orb file
				{1,2},		//number of each elements
				{12.0,1.0},	//atomic mass
				"Direct",	//coordination type
				{0.1,0.1,0.1,	//atomic coordinates
				 0.15,0.15,0.15,
				 0.05,0.05,0.05},
				{1,1,1,	//if atom can move: mbl
				 0,0,0,
				 0,0,1},
				{0.1,0.1,0.1,	//velocity: vel
				 0.1,0.1,0.1,
				 0.1,0.1,0.1}),
			UcellTestPrepare("C1H2-Cartesian",	//system-name
				"bcc",		//latname
				2,		//lmaxmax
				true,		//init_vel
				true,		//selective_dyanmics
				true,		//relax_new
				"volume",	//fixed_axes
				1.8897261254578281, //lat0
				{10.0,0.0,0.0,	//latvec
				 0.0,10.0,0.0,
				 0.0,0.0,10.0},
				{"C","H"},	//elements
				{"C.upf","H.upf"},	//upf file
				{"upf201","upf201"},	//upf types
				{"C.orb","H.orb"},	//orb file
				{1,2},		//number of each elements
				{12.0,1.0},	//atomic mass
				"Cartesian",	//coordination type
				{1,1,1,	//atomic coordinates
				 1.5,1.5,1.5,
				 0.5,0.5,0.5},
				{1,1,1,	//if atom can move: mbl
				 0,0,0,
				 0,0,1},
				{0.1,0.1,0.1,	//velocity: vel
				 0.1,0.1,0.1,
				 0.1,0.1,0.1}),
			UcellTestPrepare("C1H2-CheckDTau",	//system-name
				"bcc",		//latname
				2,		//lmaxmax
				false,		//init_vel
				false,		//selective_dyanmics
				true,		//relax_new
				"volume",	//fixed_axes
				1.8897261254578281, //lat0
				{0.1,0.1,0.1,	//latvec
				 0.15,0.15,0.15,
				 0.05,0.05,0.05},
				{"C","H"},	//elements
				{"C.upf","H.upf"},	//upf file
				{"upf201","upf201"},	//upf types
				{"C.orb","H.orb"},	//orb file
				{1,2},		//number of each elements
				{12.0,1.0},	//atomic mass
				"Direct",	//coordination type
				{1.6,2.5,3.8,	//atomic coordinates
				 -0.15,1.0,-0.15,
				 -3.05,-2.8,0.0}),
			UcellTestPrepare("C1H2-CheckTau",	//system-name
				"bcc",		//latname
				2,		//lmaxmax
				false,		//init_vel
				false,		//selective_dyanmics
				true,		//relax_new
				"volume",	//fixed_axes
				1.8897261254578281, //lat0
				{0.1,0.1,0.1,	//latvec
				 0.15,0.15,0.15,
				 0.05,0.05,0.05},
				{"C","H"},	//elements
				{"C.upf","H.upf"},	//upf file
				{"upf201","upf201"},	//upf types
				{"C.orb","H.orb"},	//orb file
				{1,2},		//number of each elements
				{12.0,1.0},	//atomic mass
				"Direct",	//coordination type
				{0.0,0.0,0.0,	//atomic coordinates
				 0.00001,0.00001,0.00001,
				 -3.05,-2.8,0.0}),
			UcellTestPrepare("C1H2-SD",	//system-name
				"bcc",		//latname
				2,		//lmaxmax
				false,		//init_vel
				false,		//selective_dyanmics
				true,		//relax_new
				"volume",	//fixed_axes
				1.8897261254578281, //lat0
				{0.1,0.1,0.1,	//latvec
				 0.15,0.15,0.15,
				 0.05,0.05,0.05},
				{"C","H"},	//elements
				{"C.upf","H.upf"},	//upf file
				{"upf201","upf201"},	//upf types
				{"C.orb","H.orb"},	//orb file
				{1,2},		//number of each elements
				{12.0,1.0},	//atomic mass
				"Direct",	//coordination type
				{0.1,0.1,0.1,	//atomic coordinates
				 0.15,0.15,0.15,
				 0.05,0.05,0.05},
				{1,1,1,	//if atom can move
				 0,0,0,
				 0,0,1},
				{0.1,0.1,0.1,	//velocity
				 0.1,0.1,0.1,
				 0.1,0.1,0.1}),
			UcellTestPrepare("C1H2-PBA",	//system-name
				"bcc",		//latname
				2,		//lmaxmax
				false,		//init_vel
				false,		//selective_dyanmics
				true,		//relax_new
				"volume",	//fixed_axes
				1.8897261254578281, //lat0
				{0.1,0.1,0.1,	//latvec
				 0.15,0.15,0.15,
				 0.05,0.05,0.05},
				{"C","H"},	//elements
				{"C.upf","H.upf"},	//upf file
				{"upf201","upf201"},	//upf types
				{"C.orb","H.orb"},	//orb file
				{1,2},		//number of each elements
				{12.0,1.0},	//atomic mass
				"Direct",	//coordination type
				{-0.1,-0.1,-0.1,	//atomic coordinates
				 1.2,1.2,1.2,
				 -3.05,-2.8,0.0})
			));
