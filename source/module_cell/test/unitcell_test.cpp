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
 *   - Constructor:
 *     - UnitCell() and ~UnitCell()
 *   - Setup:
 *     - setup(): to set latname, ntype, lmaxmax, init_vel, and lc
 *     - if_cell_can_change(): judge if any lattice vector can change
 *   - SetupWarningQuit1:
 *     - setup(): deliver warning: "there are bugs in the old implementation; 
 *         set relax_new to be 1 for fixed_volume relaxation"
 *   - SetupWarningQuit2:
 *     - setup(): deliver warning: "set relax_new to be 1 for fixed_shape relaxation"
 *   - RemakeCell
 *     - remake_cell(): rebuild cell according to its latName
 *   - RemakeCellWarnings
 *     - remake_cell(): deliver warnings when find wrong latname or cos12
 *   - JudgeParallel
 *     - judge_parallel: judge if two vectors a[3] and Vector3<double> b are parallel
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
 *   - PeriodicBoundaryAdjustment
 *     - periodic_boundary_adjustment(): move atoms inside the unitcell after relaxation
 *   - PrintCell
 *     - print_cell(ofs): print basic cell info into ofs
 *   - PrintCellCif
 *     - print_cell_cif(fn): print cif file (named as fn)
 *   - PrintSTRU
 *     - print_stru_file(): print STRU file of ABACUS
 *   - PrintTauDirect
 *   - PrintTauCartesian
 *     - print_tau(): print atomic coordinates, magmom and initial velocities
 *   - PrintUnitcellPseudo
 *     - Actually an integrated function to call UnitCell::print_cell and Atom::print_Atom
 *   - UpdateVel
 *     - update_vel(const ModuleBase::Vector3<double>* vel_in)
 *   - CalUx
 *     - cal_ux(): calculate magnetic moments of cell
 *   - ReadOrbFile
 *     - read_orb_file(): read header part of orbital file
 *   - ReadOrbFileWarning
 *     - read_orb_file(): ABACUS Cannot find the ORBITAL file
 *   - ReadAtomSpecies
 *     - read_atom_species(): a successful case
 *   - ReadAtomSpeciesWarning1
 *     - read_atom_species(): unrecongnized pseudo type.
 *   - ReadAtomSpeciesWarning2
 *     - read_atom_species(): lat0<=0.0
 *   - ReadAtomSpeciesWarning3
 *     - read_atom_species(): do not use LATTICE_PARAMETERS without explicit specification of lattice type
 *   - ReadAtomSpeciesWarning4
 *     - read_atom_species():do not use LATTICE_VECTORS along with explicit specification of lattice type
 *   - ReadAtomSpeciesWarning5
 *     - read_atom_species():latname not supported
 *   - ReadAtomSpeciesLatName
 *     - read_atom_species(): various latname
 *   - ReadAtomPositionsS1
 *     - read_atom_positions(): spin 1 case
 *   - ReadAtomPositionsS2
 *     - read_atom_positions(): spin 2 case
 *   - ReadAtomPositionsS4Noncolin
 *     - read_atom_positions(): spin 4 noncolinear case
 *   - ReadAtomPositionsS4Colin
 *     - read_atom_positions(): spin 4 colinear case
 *   - ReadAtomPositionsC
 *     - read_atom_positions(): Cartesian coordinates
 *   - ReadAtomPositionsCA
 *     - read_atom_positions(): Cartesian_angstrom coordinates
 *   - ReadAtomPositionsCACXY
 *     - read_atom_positions(): Cartesian_angstrom_center_xy coordinates
 *   - ReadAtomPositionsCACXZ
 *     - read_atom_positions(): Cartesian_angstrom_center_xz coordinates
 *   - ReadAtomPositionsCACXYZ
 *     - read_atom_positions(): Cartesian_angstrom_center_xyz coordinates
 *   - ReadAtomPositionsCAU
 *     - read_atom_positions(): Cartesian_au coordinates
 *   - ReadAtomPositionsWarning1
 *     - read_atom_positions(): unknown type of coordinates
 *   - ReadAtomPositionsWarning2
 *     - read_atom_positions(): atomic label inconsistency between ATOM_POSITIONS
 *                              and ATOM_SPECIES
 *   - ReadAtomPositionsWarning3
 *     - read_atom_positions(): warning :  atom number < 0
 *   - ReadAtomPositionsWarning4
 *     - read_atom_positions(): mismatch in atom number for atom type
 *   - ReadAtomPositionsWarning5
 *     - read_atom_positions(): no atom can move in MD!
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
	std::unique_ptr<UnitCell> ucell{new UnitCell};
	std::string output;
};

using UcellDeathTest = UcellTest;

TEST_F(UcellTest,Constructor)
{
	GlobalV::test_unitcell = 1;
	EXPECT_EQ(ucell->Coordinate,"Direct");
	EXPECT_EQ(ucell->latName,"none");
	EXPECT_DOUBLE_EQ(ucell->lat0,0.0);
	EXPECT_DOUBLE_EQ(ucell->lat0_angstrom,0.0);
	EXPECT_EQ(ucell->ntype,0);
	EXPECT_EQ(ucell->nat,0);
	EXPECT_EQ(ucell->namax,0);
	EXPECT_EQ(ucell->nwmax,0);
	EXPECT_EQ(ucell->iat2it,nullptr);
	EXPECT_EQ(ucell->iat2ia,nullptr);
	EXPECT_EQ(ucell->iwt2iat,nullptr);
	EXPECT_EQ(ucell->iwt2iw,nullptr);
	EXPECT_DOUBLE_EQ(ucell->tpiba,0.0);
	EXPECT_DOUBLE_EQ(ucell->tpiba2,0.0);
	EXPECT_DOUBLE_EQ(ucell->omega,0.0);
	EXPECT_EQ(ucell->atom_mass,nullptr);
	EXPECT_FALSE(ucell->set_atom_flag);
}

TEST_F(UcellTest,Setup)
{
	std::string latname_in = "bcc";
	int ntype_in = 1;
	int lmaxmax_in = 2;
	bool init_vel_in = false;
	std::vector<std::string> fixed_axes_in = {"None","volume","shape","a","b","c","ab","ac","bc","abc"};
	GlobalV::relax_new = true;
	for(int i=0;i<fixed_axes_in.size();++i)
	{
		ucell->setup(latname_in,ntype_in,lmaxmax_in,init_vel_in,fixed_axes_in[i]);
		EXPECT_EQ(ucell->latName,latname_in);
		EXPECT_EQ(ucell->ntype,ntype_in);
		EXPECT_EQ(ucell->lmaxmax,lmaxmax_in);
		EXPECT_EQ(ucell->init_vel,init_vel_in);
		if(fixed_axes_in[i] == "None" || fixed_axes_in[i] == "volume"
				|| fixed_axes_in[i] == "shape")
		{
			EXPECT_EQ(ucell->lc[0],1);
			EXPECT_EQ(ucell->lc[1],1);
			EXPECT_EQ(ucell->lc[2],1);
			EXPECT_TRUE(ucell->if_cell_can_change());
		}
		else if(fixed_axes_in[i] == "a")
		{
			EXPECT_EQ(ucell->lc[0],0);
			EXPECT_EQ(ucell->lc[1],1);
			EXPECT_EQ(ucell->lc[2],1);
			EXPECT_TRUE(ucell->if_cell_can_change());
		}
		else if(fixed_axes_in[i] == "b")
		{
			EXPECT_EQ(ucell->lc[0],1);
			EXPECT_EQ(ucell->lc[1],0);
			EXPECT_EQ(ucell->lc[2],1);
			EXPECT_TRUE(ucell->if_cell_can_change());
		}
		else if(fixed_axes_in[i] == "c")
		{
			EXPECT_EQ(ucell->lc[0],1);
			EXPECT_EQ(ucell->lc[1],1);
			EXPECT_EQ(ucell->lc[2],0);
			EXPECT_TRUE(ucell->if_cell_can_change());
		}
		else if(fixed_axes_in[i] == "ab")
		{
			EXPECT_EQ(ucell->lc[0],0);
			EXPECT_EQ(ucell->lc[1],0);
			EXPECT_EQ(ucell->lc[2],1);
			EXPECT_TRUE(ucell->if_cell_can_change());
		}
		else if(fixed_axes_in[i] == "ac")
		{
			EXPECT_EQ(ucell->lc[0],0);
			EXPECT_EQ(ucell->lc[1],1);
			EXPECT_EQ(ucell->lc[2],0);
			EXPECT_TRUE(ucell->if_cell_can_change());
		}
		else if(fixed_axes_in[i] == "bc")
		{
			EXPECT_EQ(ucell->lc[0],1);
			EXPECT_EQ(ucell->lc[1],0);
			EXPECT_EQ(ucell->lc[2],0);
			EXPECT_TRUE(ucell->if_cell_can_change());
		}
		else if(fixed_axes_in[i] == "abc")
		{
			EXPECT_EQ(ucell->lc[0],0);
			EXPECT_EQ(ucell->lc[1],0);
			EXPECT_EQ(ucell->lc[2],0);
			EXPECT_FALSE(ucell->if_cell_can_change());
		}
	}
}

TEST_F(UcellDeathTest,SetupWarningQuit1)
{
	std::string latname_in = "bcc";
	int ntype_in = 1;
	int lmaxmax_in = 2;
	bool init_vel_in = false;
	GlobalV::relax_new = false;
	std::string fixed_axes_in = "volume";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(ucell->setup(latname_in,ntype_in,lmaxmax_in,init_vel_in,fixed_axes_in),
			::testing::ExitedWithCode(0),"");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("there are bugs in the old implementation; set relax_new to be 1 for fixed_volume relaxation"));
}

TEST_F(UcellDeathTest,SetupWarningQuit2)
{
	std::string latname_in = "bcc";
	int ntype_in = 1;
	int lmaxmax_in = 2;
	bool init_vel_in = false;
	GlobalV::relax_new = false;
	std::string fixed_axes_in = "shape";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(ucell->setup(latname_in,ntype_in,lmaxmax_in,init_vel_in,fixed_axes_in),
			::testing::ExitedWithCode(0),"");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("set relax_new to be 1 for fixed_shape relaxation"));
}

TEST_F(UcellDeathTest, CompareAatomLabel)
{
    std::string stru_label[] =   {"Ag"  "Ag"  "Ag"     "47" "47" "47"     "Silver" "Silver" "Silver" "Ag"  "Ag" "Ag"         "Ag_empty"};
    std::string pseudo_label[] = {"Ag"  "47"  "Silver" "Ag" "47" "Silver" "Ag"     "47"     "Silver" "Ag1" "ag" "ag_locpsp"  "Ag"      };
	for (it = 0; it < 12; it++)
	{
	ucell->compare_atom_labels(stru_label[it], pseudo_label[it]);
	}
	stru_label[0] = "Fe";
	pseudo_label[0] = "O";
    testing::internal::CaptureStdout();
    EXPECT_EXIT(ucell->compare_atom_labels(stru_label[0], pseudo_label[0]),::testing::ExitedWithCode(0),"");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("atom label in STRU is " + stru_label + " mismatch with pseudo file of " +pseudo_label));
}

TEST_F(UcellTest,RemakeCell)
{
	std::vector<std::string> latname_in = {"sc","fcc","bcc","hexagonal","trigonal","st","bct","so","baco","fco","bco","sm","bacm","triclinic"};
	for(int i=0; i<latname_in.size(); ++i)
	{
		ucell->latvec.e11 = 10.0; ucell->latvec.e12 = 0.00; ucell->latvec.e13 = 0.00;
		ucell->latvec.e21 = 0.00; ucell->latvec.e22 = 10.0; ucell->latvec.e23 = 0.00;
		ucell->latvec.e31 = 0.00; ucell->latvec.e32 = 0.00; ucell->latvec.e33 = 10.0;
		ucell->latName = latname_in[i];
		ucell->remake_cell();
		if(latname_in[i]=="sc")
		{
			double celldm = std::sqrt(pow(ucell->latvec.e11,2)+pow(ucell->latvec.e12,2)+pow(ucell->latvec.e13,2));
			EXPECT_DOUBLE_EQ(ucell->latvec.e11,celldm);
		}
		else if(latname_in[i]=="fcc")
		{
			double celldm = std::sqrt(pow(ucell->latvec.e11,2)+pow(ucell->latvec.e12,2)+pow(ucell->latvec.e13,2)) / std::sqrt(2.0);
			EXPECT_DOUBLE_EQ(ucell->latvec.e11,-celldm);
			EXPECT_DOUBLE_EQ(ucell->latvec.e12,0.0);
			EXPECT_DOUBLE_EQ(ucell->latvec.e13,celldm);
			EXPECT_DOUBLE_EQ(ucell->latvec.e21,0.0);
			EXPECT_DOUBLE_EQ(ucell->latvec.e22,celldm);
			EXPECT_DOUBLE_EQ(ucell->latvec.e23,celldm);
			EXPECT_DOUBLE_EQ(ucell->latvec.e31,-celldm);
			EXPECT_DOUBLE_EQ(ucell->latvec.e32,celldm);
			EXPECT_DOUBLE_EQ(ucell->latvec.e33,0.0);
		}
		else if(latname_in[i]=="bcc")
		{
			double celldm = std::sqrt(pow(ucell->latvec.e11,2)+pow(ucell->latvec.e12,2)+pow(ucell->latvec.e13,2)) / std::sqrt(3.0);
			EXPECT_DOUBLE_EQ(ucell->latvec.e11,celldm);
			EXPECT_DOUBLE_EQ(ucell->latvec.e12,celldm);
			EXPECT_DOUBLE_EQ(ucell->latvec.e13,celldm);
			EXPECT_DOUBLE_EQ(ucell->latvec.e21,-celldm);
			EXPECT_DOUBLE_EQ(ucell->latvec.e22,celldm);
			EXPECT_DOUBLE_EQ(ucell->latvec.e23,celldm);
			EXPECT_DOUBLE_EQ(ucell->latvec.e31,-celldm);
			EXPECT_DOUBLE_EQ(ucell->latvec.e32,-celldm);
			EXPECT_DOUBLE_EQ(ucell->latvec.e33,celldm);
		}
		else if(latname_in[i]=="hexagonal")
		{
			double celldm1 = std::sqrt(pow(ucell->latvec.e11,2)+pow(ucell->latvec.e12,2)+pow(ucell->latvec.e13,2));
			double celldm3 = std::sqrt(pow(ucell->latvec.e31,2)+pow(ucell->latvec.e32,2)+pow(ucell->latvec.e33,2));
			double mathfoo = sqrt(3.0)/2.0;
			EXPECT_DOUBLE_EQ(ucell->latvec.e11,celldm1);
			EXPECT_DOUBLE_EQ(ucell->latvec.e12,0.0);
			EXPECT_DOUBLE_EQ(ucell->latvec.e13,0.0);
			EXPECT_DOUBLE_EQ(ucell->latvec.e21,-0.5*celldm1);
			EXPECT_DOUBLE_EQ(ucell->latvec.e22,celldm1*mathfoo);
			EXPECT_DOUBLE_EQ(ucell->latvec.e23,0.0);
			EXPECT_DOUBLE_EQ(ucell->latvec.e31,0.0);
			EXPECT_DOUBLE_EQ(ucell->latvec.e32,0.0);
			EXPECT_DOUBLE_EQ(ucell->latvec.e33,celldm3);
		}
		else if(latname_in[i]=="trigonal")
		{
			double a1 = std::sqrt(pow(ucell->latvec.e11,2)+pow(ucell->latvec.e12,2)+pow(ucell->latvec.e13,2));
			double a2 = std::sqrt(pow(ucell->latvec.e21,2)+pow(ucell->latvec.e22,2)+pow(ucell->latvec.e23,2));
			double a1da2 = (ucell->latvec.e11*ucell->latvec.e21 + ucell->latvec.e12*ucell->latvec.e22 + ucell->latvec.e13*ucell->latvec.e23);
			double cosgamma = a1da2 / (a1 * a2);
			double tx = std::sqrt((1.0 - cosgamma)/2.0);
			double ty = std::sqrt((1.0 - cosgamma)/6.0);
			double tz = std::sqrt((1.0 + 2.0*cosgamma)/3.0);
			EXPECT_DOUBLE_EQ(ucell->latvec.e11, a1*tx); EXPECT_DOUBLE_EQ(ucell->latvec.e12,-a1*ty); EXPECT_DOUBLE_EQ(ucell->latvec.e13,a1*tz);
			EXPECT_DOUBLE_EQ(ucell->latvec.e21, 0.0); EXPECT_DOUBLE_EQ(ucell->latvec.e22,2.0*a1*ty); EXPECT_DOUBLE_EQ(ucell->latvec.e23,a1*tz);
			EXPECT_DOUBLE_EQ(ucell->latvec.e31, -a1*tx); EXPECT_DOUBLE_EQ(ucell->latvec.e32,-a1*ty); EXPECT_DOUBLE_EQ(ucell->latvec.e33,a1*tz);
		}
		else if(latname_in[i]=="st")
		{
			double a1 = std::sqrt(pow(ucell->latvec.e11,2)+pow(ucell->latvec.e12,2)+pow(ucell->latvec.e13,2));
			double a3 = std::sqrt(pow(ucell->latvec.e31,2)+pow(ucell->latvec.e32,2)+pow(ucell->latvec.e33,2));
			EXPECT_DOUBLE_EQ(ucell->latvec.e11, a1); EXPECT_DOUBLE_EQ(ucell->latvec.e12,0.0); EXPECT_DOUBLE_EQ(ucell->latvec.e13,0.0);
			EXPECT_DOUBLE_EQ(ucell->latvec.e21, 0.0); EXPECT_DOUBLE_EQ(ucell->latvec.e22,a1); EXPECT_DOUBLE_EQ(ucell->latvec.e23,0.0);
			EXPECT_DOUBLE_EQ(ucell->latvec.e31, 0.0); EXPECT_DOUBLE_EQ(ucell->latvec.e32,0.0); EXPECT_DOUBLE_EQ(ucell->latvec.e33,a3);
		}
		else if(latname_in[i]=="bct")
		{
			double d1 = std::abs(ucell->latvec.e11);
			double d2 = std::abs(ucell->latvec.e13);
			EXPECT_DOUBLE_EQ(ucell->latvec.e11, d1); EXPECT_DOUBLE_EQ(ucell->latvec.e12,-d1); EXPECT_DOUBLE_EQ(ucell->latvec.e13,d2);
			EXPECT_DOUBLE_EQ(ucell->latvec.e21, d1); EXPECT_DOUBLE_EQ(ucell->latvec.e22, d1); EXPECT_DOUBLE_EQ(ucell->latvec.e23,d2);
			EXPECT_DOUBLE_EQ(ucell->latvec.e31,-d1); EXPECT_DOUBLE_EQ(ucell->latvec.e32,-d1); EXPECT_DOUBLE_EQ(ucell->latvec.e33,d2);
		}
		else if(latname_in[i]=="so")
		{
			double a1 = std::sqrt(pow(ucell->latvec.e11,2)+pow(ucell->latvec.e12,2)+pow(ucell->latvec.e13,2));
			double a2 = std::sqrt(pow(ucell->latvec.e21,2)+pow(ucell->latvec.e22,2)+pow(ucell->latvec.e23,2));
			double a3 = std::sqrt(pow(ucell->latvec.e31,2)+pow(ucell->latvec.e32,2)+pow(ucell->latvec.e33,2));
			EXPECT_DOUBLE_EQ(ucell->latvec.e11, a1); EXPECT_DOUBLE_EQ(ucell->latvec.e12,0.0); EXPECT_DOUBLE_EQ(ucell->latvec.e13,0.0);
			EXPECT_DOUBLE_EQ(ucell->latvec.e21,0.0); EXPECT_DOUBLE_EQ(ucell->latvec.e22,a2); EXPECT_DOUBLE_EQ(ucell->latvec.e23,0.0);
			EXPECT_DOUBLE_EQ(ucell->latvec.e31,0.0); EXPECT_DOUBLE_EQ(ucell->latvec.e32,0.0); EXPECT_DOUBLE_EQ(ucell->latvec.e33,a3);
		}
		else if(latname_in[i]=="baco")
		{
			double d1 = std::abs(ucell->latvec.e11);
			double d2 = std::abs(ucell->latvec.e22);
			double d3 = std::abs(ucell->latvec.e33);
			EXPECT_DOUBLE_EQ(ucell->latvec.e11, d1); EXPECT_DOUBLE_EQ(ucell->latvec.e12,d2); EXPECT_DOUBLE_EQ(ucell->latvec.e13,0.0);
			EXPECT_DOUBLE_EQ(ucell->latvec.e21,-d1); EXPECT_DOUBLE_EQ(ucell->latvec.e22,d2); EXPECT_DOUBLE_EQ(ucell->latvec.e23,0.0);
			EXPECT_DOUBLE_EQ(ucell->latvec.e31,0.0); EXPECT_DOUBLE_EQ(ucell->latvec.e32,0.0); EXPECT_DOUBLE_EQ(ucell->latvec.e33,d3);
		}
		else if(latname_in[i]=="fco")
		{
			double d1 = std::abs(ucell->latvec.e11);
			double d2 = std::abs(ucell->latvec.e22);
			double d3 = std::abs(ucell->latvec.e33);
			EXPECT_DOUBLE_EQ(ucell->latvec.e11, d1); EXPECT_DOUBLE_EQ(ucell->latvec.e12,0.0); EXPECT_DOUBLE_EQ(ucell->latvec.e13,d3);
			EXPECT_DOUBLE_EQ(ucell->latvec.e21, d1); EXPECT_DOUBLE_EQ(ucell->latvec.e22,d2); EXPECT_DOUBLE_EQ(ucell->latvec.e23,0.0);
			EXPECT_DOUBLE_EQ(ucell->latvec.e31,0.0); EXPECT_DOUBLE_EQ(ucell->latvec.e32,d2); EXPECT_DOUBLE_EQ(ucell->latvec.e33,d3);
		}
		else if(latname_in[i]=="bco")
		{
			double d1 = std::abs(ucell->latvec.e11);
			double d2 = std::abs(ucell->latvec.e22);
			double d3 = std::abs(ucell->latvec.e33);
			EXPECT_DOUBLE_EQ(ucell->latvec.e11, d1); EXPECT_DOUBLE_EQ(ucell->latvec.e12,d2); EXPECT_DOUBLE_EQ(ucell->latvec.e13,d3);
			EXPECT_DOUBLE_EQ(ucell->latvec.e21,-d1); EXPECT_DOUBLE_EQ(ucell->latvec.e22,d2); EXPECT_DOUBLE_EQ(ucell->latvec.e23,d3);
			EXPECT_DOUBLE_EQ(ucell->latvec.e31,-d1); EXPECT_DOUBLE_EQ(ucell->latvec.e32,-d2); EXPECT_DOUBLE_EQ(ucell->latvec.e33,d3);
		}
		else if(latname_in[i]=="sm")
		{
			double a1 = std::sqrt(pow(ucell->latvec.e11,2)+pow(ucell->latvec.e12,2)+pow(ucell->latvec.e13,2));
			double a2 = std::sqrt(pow(ucell->latvec.e21,2)+pow(ucell->latvec.e22,2)+pow(ucell->latvec.e23,2));
			double a3 = std::sqrt(pow(ucell->latvec.e31,2)+pow(ucell->latvec.e32,2)+pow(ucell->latvec.e33,2));
			double a1da2 = (ucell->latvec.e11*ucell->latvec.e21 + ucell->latvec.e12*ucell->latvec.e22 + ucell->latvec.e13*ucell->latvec.e23);
			double cosgamma = a1da2 / (a1 * a2);
			double d1 = a2 * cosgamma;
			double d2 = a2 * std::sqrt(1.0 - cosgamma * cosgamma);
			EXPECT_DOUBLE_EQ(ucell->latvec.e11, a1); EXPECT_DOUBLE_EQ(ucell->latvec.e12,0.0); EXPECT_DOUBLE_EQ(ucell->latvec.e13,0.0);
			EXPECT_DOUBLE_EQ(ucell->latvec.e21, d1); EXPECT_DOUBLE_EQ(ucell->latvec.e22,d2); EXPECT_DOUBLE_EQ(ucell->latvec.e23,0.0);
			EXPECT_DOUBLE_EQ(ucell->latvec.e31, 0.0); EXPECT_DOUBLE_EQ(ucell->latvec.e32,0.0); EXPECT_DOUBLE_EQ(ucell->latvec.e33,a3);
		}
		else if(latname_in[i]=="bacm")
		{
			double d1 = std::abs(ucell->latvec.e11);
			double a2 = std::sqrt(pow(ucell->latvec.e21,2)+pow(ucell->latvec.e22,2)+pow(ucell->latvec.e23,2));
			double d3 = std::abs(ucell->latvec.e13);
			double cosgamma = ucell->latvec.e21/a2;
			double f1 = a2 * cosgamma;
			double f2 = a2 * std::sqrt(1.0 - cosgamma*cosgamma);
			EXPECT_DOUBLE_EQ(ucell->latvec.e11, d1); EXPECT_DOUBLE_EQ(ucell->latvec.e12,0.0); EXPECT_DOUBLE_EQ(ucell->latvec.e13,-d3);
			EXPECT_DOUBLE_EQ(ucell->latvec.e21, f1); EXPECT_DOUBLE_EQ(ucell->latvec.e22,f2); EXPECT_DOUBLE_EQ(ucell->latvec.e23,0.0);
			EXPECT_DOUBLE_EQ(ucell->latvec.e31, d1); EXPECT_DOUBLE_EQ(ucell->latvec.e32,0.0); EXPECT_DOUBLE_EQ(ucell->latvec.e33,d3);
		}
		else if(latname_in[i]=="triclinic")
		{
			double a1 = std::sqrt(pow(ucell->latvec.e11,2)+pow(ucell->latvec.e12,2)+pow(ucell->latvec.e13,2));
			double a2 = std::sqrt(pow(ucell->latvec.e21,2)+pow(ucell->latvec.e22,2)+pow(ucell->latvec.e23,2));
			double a3 = std::sqrt(pow(ucell->latvec.e31,2)+pow(ucell->latvec.e32,2)+pow(ucell->latvec.e33,2));
			double a1da2 = (ucell->latvec.e11*ucell->latvec.e21 + ucell->latvec.e12*ucell->latvec.e22 + ucell->latvec.e13*ucell->latvec.e23);
			double a1da3 = (ucell->latvec.e11*ucell->latvec.e31 + ucell->latvec.e12*ucell->latvec.e32 + ucell->latvec.e13*ucell->latvec.e33);
			double a2da3 = (ucell->latvec.e21*ucell->latvec.e31 + ucell->latvec.e22*ucell->latvec.e32 + ucell->latvec.e23*ucell->latvec.e33);
			double cosgamma = a1da2 / a1 / a2;
			double singamma = std::sqrt(1.0 - cosgamma*cosgamma);
			double cosbeta = a1da3 / a1 / a3;
			double cosalpha = a2da3 / a2 / a3;
			double d1 = std::sqrt(1.0 + 2.0*cosgamma*cosbeta*cosalpha - cosgamma*cosgamma - cosbeta*cosbeta - cosalpha*cosalpha)/singamma;
			EXPECT_DOUBLE_EQ(ucell->latvec.e11, a1); EXPECT_DOUBLE_EQ(ucell->latvec.e12,0.0); EXPECT_DOUBLE_EQ(ucell->latvec.e13,0.0);
			EXPECT_DOUBLE_EQ(ucell->latvec.e21, a2*cosgamma); EXPECT_DOUBLE_EQ(ucell->latvec.e22,a2*singamma); EXPECT_DOUBLE_EQ(ucell->latvec.e23,0.0);
			EXPECT_DOUBLE_EQ(ucell->latvec.e31, a3*cosbeta); EXPECT_DOUBLE_EQ(ucell->latvec.e32,a3*(cosalpha - cosbeta*cosgamma)/singamma); 
			EXPECT_DOUBLE_EQ(ucell->latvec.e33, a3*d1);
		}
	}
}

TEST_F(UcellDeathTest,RemakeCellWarnings)
{
	std::vector<std::string> latname_in = {"none","trigonal","bacm","triclinic","arbitrary"};
	for(int i=0; i<latname_in.size(); ++i)
	{
		ucell->latvec.e11 = 10.0; ucell->latvec.e12 = 0.00; ucell->latvec.e13 = 0.00;
		ucell->latvec.e21 = 10.0; ucell->latvec.e22 = 0.00; ucell->latvec.e23 = 0.00;
		ucell->latvec.e31 = 0.00; ucell->latvec.e32 = 0.00; ucell->latvec.e33 = 10.0;
		ucell->latName = latname_in[i];
		testing::internal::CaptureStdout();
		EXPECT_EXIT(ucell->remake_cell(),::testing::ExitedWithCode(0),"");
		std::string output = testing::internal::GetCapturedStdout();
		if(latname_in[i]=="none")
		{
			EXPECT_THAT(output,testing::HasSubstr("to use fixed_ibrav, latname must be provided"));
		}
		else if(latname_in[i]=="trigonal" || latname_in[i]=="bacm" || latname_in[i] == "triclinic")
		{
			EXPECT_THAT(output,testing::HasSubstr("wrong cos12!"));
		}
		else
		{
			EXPECT_THAT(output,testing::HasSubstr("latname not supported!"));
		}
	}
}


TEST_F(UcellTest,JudgeParallel)
{
	ModuleBase::Vector3<double> b(1.0,1.0,1.0);
	double a[3] = {1.0,1.0,1.0};
	EXPECT_TRUE(ucell->judge_parallel(a,b));
}

TEST_F(UcellTest,Index)
{
	UcellTestPrepare utp = UcellTestLib["C1H2-Index"];
	GlobalV::relax_new = utp.relax_new;
	ucell = utp.SetUcellInfo();
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

TEST_F(UcellTest,CheckDTau)
{
	UcellTestPrepare utp = UcellTestLib["C1H2-CheckDTau"];
	GlobalV::relax_new = utp.relax_new;
	ucell = utp.SetUcellInfo();	
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

TEST_F(UcellTest,CheckTau)
{
	UcellTestPrepare utp=UcellTestLib["C1H2-CheckTau"];
	GlobalV::relax_new = utp.relax_new;
	ucell = utp.SetUcellInfo();
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


TEST_F(UcellTest,SelectiveDynamics)
{
	UcellTestPrepare utp=UcellTestLib["C1H2-SD"];
	GlobalV::relax_new = utp.relax_new;
	ucell = utp.SetUcellInfo();
	EXPECT_TRUE(ucell->if_atoms_can_move());
}

TEST_F(UcellDeathTest,PeriodicBoundaryAdjustment1)
{
	UcellTestPrepare utp = UcellTestLib["C1H2-PBA"];
	GlobalV::relax_new = utp.relax_new;
	ucell = utp.SetUcellInfo();
	testing::internal::CaptureStdout();
	EXPECT_EXIT(ucell->periodic_boundary_adjustment(),::testing::ExitedWithCode(0),"");
	std::string output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("the movement of atom is larger than the length of cell"));
}

TEST_F(UcellTest,PeriodicBoundaryAdjustment2)
{
	UcellTestPrepare utp = UcellTestLib["C1H2-Index"];
	GlobalV::relax_new = utp.relax_new;
	ucell = utp.SetUcellInfo();
	EXPECT_NO_THROW(ucell->periodic_boundary_adjustment());
}

TEST_F(UcellTest,PrintCell)
{
	UcellTestPrepare utp = UcellTestLib["C1H2-Index"];
	GlobalV::relax_new = utp.relax_new;
	ucell = utp.SetUcellInfo();
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

TEST_F(UcellTest,PrintUnitcellPseudo)
{
	UcellTestPrepare utp = UcellTestLib["C1H2-Index"];
	GlobalV::relax_new = utp.relax_new;
	ucell = utp.SetUcellInfo();
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

TEST_F(UcellTest,PrintCellCif)
{
	UcellTestPrepare utp = UcellTestLib["C1H2-Index"];
	GlobalV::relax_new = utp.relax_new;
	ucell = utp.SetUcellInfo();
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

TEST_F(UcellTest,PrintSTRU)
{
	UcellTestPrepare utp = UcellTestLib["C1H2-Index"];
	GlobalV::relax_new = utp.relax_new;
	ucell = utp.SetUcellInfo();
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
    	EXPECT_THAT(str, testing::HasSubstr("1.5000000000     1.5000000000     1.5000000000 m  0  0  0 v     0.1000000000     0.1000000000     0.1000000000"));
    	EXPECT_THAT(str, testing::HasSubstr("0.5000000000     0.5000000000     0.5000000000 m  0  0  1 v     0.1000000000     0.1000000000     0.1000000000"));
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
    	EXPECT_THAT(str, testing::HasSubstr("0.1500000000     0.1500000000     0.1500000000 m  0  0  0 v     0.1000000000     0.1000000000     0.1000000000"));
    	EXPECT_THAT(str, testing::HasSubstr("0.0500000000     0.0500000000     0.0500000000 m  0  0  1 v     0.1000000000     0.1000000000     0.1000000000"));
	ifs.close();
	remove("C1H2_STRU");
}

TEST_F(UcellTest,PrintTauDirect)
{
	UcellTestPrepare utp = UcellTestLib["C1H2-Index"];
	GlobalV::relax_new = utp.relax_new;
	ucell = utp.SetUcellInfo();
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

TEST_F(UcellTest,PrintTauCartesian)
{
	UcellTestPrepare utp = UcellTestLib["C1H2-Cartesian"];
	GlobalV::relax_new = utp.relax_new;
	ucell = utp.SetUcellInfo();
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

TEST_F(UcellTest,UpdateVel)
{
	UcellTestPrepare utp = UcellTestLib["C1H2-Index"];
	GlobalV::relax_new = utp.relax_new;
	ucell = utp.SetUcellInfo();
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

TEST_F(UcellTest,CalUx1)
{
	UcellTestPrepare utp = UcellTestLib["C1H2-Read"];
	GlobalV::relax_new = utp.relax_new;
	ucell = utp.SetUcellInfo();
	ucell->atoms[0].m_loc_[0].set(0,-1,0);
	ucell->atoms[1].m_loc_[0].set(1,1,1);
	ucell->atoms[1].m_loc_[1].set(0,0,0);
	ucell->cal_ux();
	EXPECT_FALSE(ucell->magnet.lsign_);
	EXPECT_DOUBLE_EQ(ucell->magnet.ux_[0],0);
	EXPECT_DOUBLE_EQ(ucell->magnet.ux_[1],-1);
	EXPECT_DOUBLE_EQ(ucell->magnet.ux_[2],0);
}

TEST_F(UcellTest,CalUx2)
{
	UcellTestPrepare utp = UcellTestLib["C1H2-Read"];
	GlobalV::relax_new = utp.relax_new;
	ucell = utp.SetUcellInfo();
	ucell->atoms[0].m_loc_[0].set(0,0,0);
	ucell->atoms[1].m_loc_[0].set(1,1,1);
	ucell->atoms[1].m_loc_[1].set(0,0,0);
	//(0,0,0) is also parallel to (1,1,1)
	ucell->cal_ux();
	EXPECT_TRUE(ucell->magnet.lsign_);
	EXPECT_NEAR(ucell->magnet.ux_[0],0.57735,1e-5);
	EXPECT_NEAR(ucell->magnet.ux_[1],0.57735,1e-5);
	EXPECT_NEAR(ucell->magnet.ux_[2],0.57735,1e-5);
}


#ifdef __LCAO
TEST_F(UcellTest,ReadOrbFile)
{
	UcellTestPrepare utp = UcellTestLib["C1H2-Read"];
	GlobalV::relax_new = utp.relax_new;
	ucell = utp.SetUcellInfo();
	std::string orb_file = "./support/C.orb";
	std::ofstream ofs_running;
	ofs_running.open("tmp_readorbfile");
	ucell->read_orb_file(0,orb_file,ofs_running,&(ucell->atoms[0]));
	ofs_running.close();
	EXPECT_EQ(ucell->atoms[0].nw,25);
	remove("tmp_readorbfile");
}

TEST_F(UcellDeathTest,ReadOrbFileWarning)
{
	UcellTestPrepare utp = UcellTestLib["C1H2-Read"];
	GlobalV::relax_new = utp.relax_new;
	ucell = utp.SetUcellInfo();
	std::string orb_file = "./support/CC.orb";
	std::ofstream ofs_running;
	ofs_running.open("tmp_readorbfile");
	testing::internal::CaptureStdout();
	EXPECT_EXIT(ucell->read_orb_file(0,orb_file,ofs_running,&(ucell->atoms[0])),
			::testing::ExitedWithCode(0),"");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("ABACUS Cannot find the ORBITAL file"));
	ofs_running.close();
	remove("tmp_readorbfile");
}

TEST_F(UcellTest,ReadAtomSpecies)
{
	std::string fn = "./support/STRU_MgO";
	std::ifstream ifa(fn.c_str());
	std::ofstream ofs_running;
	ofs_running.open("read_atom_species.tmp");
	ucell->ntype = 2;
	ucell->atoms = new Atom[ucell->ntype];
	ucell->set_atom_flag = true;
	GlobalV::test_pseudo_cell = 2;
	GlobalV::BASIS_TYPE = "lcao";
	GlobalV::deepks_setorb = true;
	EXPECT_NO_THROW(ucell->read_atom_species(ifa,ofs_running));
	EXPECT_DOUBLE_EQ(ucell->latvec.e11,4.27957);
	EXPECT_DOUBLE_EQ(ucell->latvec.e22,4.27957);
	EXPECT_DOUBLE_EQ(ucell->latvec.e33,4.27957);
	ofs_running.close();
	ifa.close();
	remove("read_atom_species.tmp");
}

TEST_F(UcellDeathTest,ReadAtomSpeciesWarning1)
{
	std::string fn = "./support/STRU_MgO_Warning1";
	std::ifstream ifa(fn.c_str());
	std::ofstream ofs_running;
	ofs_running.open("read_atom_species.tmp");
	ucell->ntype = 2;
	ucell->atoms = new Atom[ucell->ntype];
	ucell->set_atom_flag = true;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(ucell->read_atom_species(ifa,ofs_running),
			::testing::ExitedWithCode(0),"");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("unrecongnized pseudo type."));
	ofs_running.close();
	ifa.close();
	remove("read_atom_species.tmp");
}

TEST_F(UcellDeathTest,ReadAtomSpeciesWarning2)
{
	std::string fn = "./support/STRU_MgO_Warning2";
	std::ifstream ifa(fn.c_str());
	std::ofstream ofs_running;
	ofs_running.open("read_atom_species.tmp");
	ucell->ntype = 2;
	ucell->atoms = new Atom[ucell->ntype];
	ucell->set_atom_flag = true;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(ucell->read_atom_species(ifa,ofs_running),
			::testing::ExitedWithCode(0),"");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("lat0<=0.0"));
	ofs_running.close();
	ifa.close();
	remove("read_atom_species.tmp");
}

TEST_F(UcellDeathTest,ReadAtomSpeciesWarning3)
{
	std::string fn = "./support/STRU_MgO_Warning3";
	std::ifstream ifa(fn.c_str());
	std::ofstream ofs_running;
	ofs_running.open("read_atom_species.tmp");
	ucell->ntype = 2;
	ucell->atoms = new Atom[ucell->ntype];
	ucell->set_atom_flag = true;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(ucell->read_atom_species(ifa,ofs_running),
			::testing::ExitedWithCode(0),"");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("do not use LATTICE_PARAMETERS without explicit specification of lattice type"));
	ofs_running.close();
	ifa.close();
	remove("read_atom_species.tmp");
}

TEST_F(UcellDeathTest,ReadAtomSpeciesWarning4)
{
	std::string fn = "./support/STRU_MgO_Warning4";
	std::ifstream ifa(fn.c_str());
	std::ofstream ofs_running;
	ofs_running.open("read_atom_species.tmp");
	ucell->ntype = 2;
	ucell->atoms = new Atom[ucell->ntype];
	ucell->set_atom_flag = true;
	ucell->latName = "bcc";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(ucell->read_atom_species(ifa,ofs_running),
			::testing::ExitedWithCode(0),"");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("do not use LATTICE_VECTORS along with explicit specification of lattice type"));
	ofs_running.close();
	ifa.close();
	remove("read_atom_species.tmp");
}


TEST_F(UcellTest,ReadAtomSpeciesLatName)
{
	ucell->ntype = 2;
	ucell->atoms = new Atom[ucell->ntype];
	ucell->set_atom_flag = true;
	std::vector<std::string> latName_in = {"sc","fcc","bcc","hexagonal","trigonal","st","bct","so","baco","fco","bco","sm","bacm","triclinic"};
	for(int i=0;i<latName_in.size();++i)
	{
		std::string fn = "./support/STRU_MgO_LatName";
		std::ifstream ifa(fn.c_str());
		std::ofstream ofs_running;
		ofs_running.open("read_atom_species.tmp");
		ucell->latName = latName_in[i];
		EXPECT_NO_THROW(ucell->read_atom_species(ifa,ofs_running));
		if(ucell->latName == "sc"){
			EXPECT_DOUBLE_EQ(ucell->latvec.e11,1.0);
			EXPECT_DOUBLE_EQ(ucell->latvec.e22,1.0);
			EXPECT_DOUBLE_EQ(ucell->latvec.e33,1.0);
		}
		ofs_running.close();
		ifa.close();
		remove("read_atom_species.tmp");
	}
}

TEST_F(UcellDeathTest,ReadAtomSpeciesWarning5)
{
	std::string fn = "./support/STRU_MgO_LatName";
	std::ifstream ifa(fn.c_str());
	std::ofstream ofs_running;
	ofs_running.open("read_atom_species.tmp");
	ucell->ntype = 2;
	ucell->atoms = new Atom[ucell->ntype];
	ucell->set_atom_flag = true;
	ucell->latName = "arbitrary";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(ucell->read_atom_species(ifa,ofs_running),
			::testing::ExitedWithCode(0),"");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("latname not supported"));
	ofs_running.close();
	ifa.close();
	remove("read_atom_species.tmp");
}

TEST_F(UcellTest,ReadAtomPositionsS1)
{
	std::string fn = "./support/STRU_MgO";
	std::ifstream ifa(fn.c_str());
	std::ofstream ofs_running;
	std::ofstream ofs_warning;
	ofs_running.open("read_atom_positions.tmp");
	ofs_warning.open("read_atom_positions.warn");
	//mandatory preliminaries
	ucell->ntype = 2;
	ucell->atoms = new Atom[ucell->ntype];
	ucell->set_atom_flag = true;
	GlobalV::test_pseudo_cell = 2;
	GlobalV::BASIS_TYPE = "lcao";
	GlobalV::deepks_setorb = true;
	GlobalV::NSPIN = 1;
	EXPECT_NO_THROW(ucell->read_atom_species(ifa,ofs_running));
	EXPECT_DOUBLE_EQ(ucell->latvec.e11,4.27957);
	EXPECT_DOUBLE_EQ(ucell->latvec.e22,4.27957);
	EXPECT_DOUBLE_EQ(ucell->latvec.e33,4.27957);
	//mandatory preliminaries
	delete[] ucell->magnet.start_magnetization;
	ucell->magnet.start_magnetization = new double[ucell->ntype];
	ucell->read_atom_positions(ifa,ofs_running,ofs_warning);
	ofs_running.close();
	ofs_warning.close();
	ifa.close();
	remove("read_atom_positions.tmp");
	remove("read_atom_positions.warn");
}

TEST_F(UcellTest,ReadAtomPositionsS2)
{
	std::string fn = "./support/STRU_MgO";
	std::ifstream ifa(fn.c_str());
	std::ofstream ofs_running;
	std::ofstream ofs_warning;
	ofs_running.open("read_atom_positions.tmp");
	ofs_warning.open("read_atom_positions.warn");
	//mandatory preliminaries
	ucell->ntype = 2;
	ucell->atoms = new Atom[ucell->ntype];
	ucell->set_atom_flag = true;
	GlobalV::test_pseudo_cell = 2;
	GlobalV::BASIS_TYPE = "lcao";
	GlobalV::deepks_setorb = true;
	GlobalV::NSPIN = 2;
	EXPECT_NO_THROW(ucell->read_atom_species(ifa,ofs_running));
	EXPECT_DOUBLE_EQ(ucell->latvec.e11,4.27957);
	EXPECT_DOUBLE_EQ(ucell->latvec.e22,4.27957);
	EXPECT_DOUBLE_EQ(ucell->latvec.e33,4.27957);
	//mandatory preliminaries
	delete[] ucell->magnet.start_magnetization;
	ucell->magnet.start_magnetization = new double[ucell->ntype];
	ucell->read_atom_positions(ifa,ofs_running,ofs_warning);
	ofs_running.close();
	ofs_warning.close();
	ifa.close();
	remove("read_atom_positions.tmp");
	remove("read_atom_positions.warn");
}

TEST_F(UcellTest,ReadAtomPositionsS4Noncolin)
{
	std::string fn = "./support/STRU_MgO";
	std::ifstream ifa(fn.c_str());
	std::ofstream ofs_running;
	std::ofstream ofs_warning;
	ofs_running.open("read_atom_positions.tmp");
	ofs_warning.open("read_atom_positions.warn");
	//mandatory preliminaries
	ucell->ntype = 2;
	ucell->atoms = new Atom[ucell->ntype];
	ucell->set_atom_flag = true;
	GlobalV::test_pseudo_cell = 2;
	GlobalV::BASIS_TYPE = "lcao";
	GlobalV::deepks_setorb = true;
	GlobalV::NSPIN = 4;
	GlobalV::NONCOLIN = true;
	EXPECT_NO_THROW(ucell->read_atom_species(ifa,ofs_running));
	EXPECT_DOUBLE_EQ(ucell->latvec.e11,4.27957);
	EXPECT_DOUBLE_EQ(ucell->latvec.e22,4.27957);
	EXPECT_DOUBLE_EQ(ucell->latvec.e33,4.27957);
	//mandatory preliminaries
	delete[] ucell->magnet.start_magnetization;
	ucell->magnet.start_magnetization = new double[ucell->ntype];
	ucell->read_atom_positions(ifa,ofs_running,ofs_warning);
	ofs_running.close();
	ofs_warning.close();
	ifa.close();
	remove("read_atom_positions.tmp");
	remove("read_atom_positions.warn");
}

TEST_F(UcellTest,ReadAtomPositionsS4Colin)
{
	std::string fn = "./support/STRU_MgO";
	std::ifstream ifa(fn.c_str());
	std::ofstream ofs_running;
	std::ofstream ofs_warning;
	ofs_running.open("read_atom_positions.tmp");
	ofs_warning.open("read_atom_positions.warn");
	//mandatory preliminaries
	ucell->ntype = 2;
	ucell->atoms = new Atom[ucell->ntype];
	ucell->set_atom_flag = true;
	GlobalV::test_pseudo_cell = 2;
	GlobalV::BASIS_TYPE = "lcao";
	GlobalV::deepks_setorb = true;
	GlobalV::NSPIN = 4;
	GlobalV::NONCOLIN = false;
	EXPECT_NO_THROW(ucell->read_atom_species(ifa,ofs_running));
	EXPECT_DOUBLE_EQ(ucell->latvec.e11,4.27957);
	EXPECT_DOUBLE_EQ(ucell->latvec.e22,4.27957);
	EXPECT_DOUBLE_EQ(ucell->latvec.e33,4.27957);
	//mandatory preliminaries
	delete[] ucell->magnet.start_magnetization;
	ucell->magnet.start_magnetization = new double[ucell->ntype];
	ucell->read_atom_positions(ifa,ofs_running,ofs_warning);
	ofs_running.close();
	ofs_warning.close();
	ifa.close();
	remove("read_atom_positions.tmp");
	remove("read_atom_positions.warn");
}

TEST_F(UcellTest,ReadAtomPositionsC)
{
	std::string fn = "./support/STRU_MgO_c";
	std::ifstream ifa(fn.c_str());
	std::ofstream ofs_running;
	std::ofstream ofs_warning;
	ofs_running.open("read_atom_positions.tmp");
	ofs_warning.open("read_atom_positions.warn");
	//mandatory preliminaries
	ucell->ntype = 2;
	ucell->atoms = new Atom[ucell->ntype];
	ucell->set_atom_flag = true;
	GlobalV::test_pseudo_cell = 2;
	GlobalV::BASIS_TYPE = "lcao";
	GlobalV::deepks_setorb = true;
	GlobalV::NSPIN = 1;
	EXPECT_NO_THROW(ucell->read_atom_species(ifa,ofs_running));
	EXPECT_DOUBLE_EQ(ucell->latvec.e11,4.27957);
	EXPECT_DOUBLE_EQ(ucell->latvec.e22,4.27957);
	EXPECT_DOUBLE_EQ(ucell->latvec.e33,4.27957);
	//mandatory preliminaries
	delete[] ucell->magnet.start_magnetization;
	ucell->magnet.start_magnetization = new double[ucell->ntype];
	ucell->read_atom_positions(ifa,ofs_running,ofs_warning);
	ofs_running.close();
	ofs_warning.close();
	ifa.close();
	remove("read_atom_positions.tmp");
	remove("read_atom_positions.warn");
}

TEST_F(UcellTest,ReadAtomPositionsCA)
{
	std::string fn = "./support/STRU_MgO_ca";
	std::ifstream ifa(fn.c_str());
	std::ofstream ofs_running;
	std::ofstream ofs_warning;
	ofs_running.open("read_atom_positions.tmp");
	ofs_warning.open("read_atom_positions.warn");
	//mandatory preliminaries
	ucell->ntype = 2;
	ucell->atoms = new Atom[ucell->ntype];
	ucell->set_atom_flag = true;
	GlobalV::test_pseudo_cell = 2;
	GlobalV::BASIS_TYPE = "lcao";
	GlobalV::deepks_setorb = true;
	GlobalV::NSPIN = 1;
	EXPECT_NO_THROW(ucell->read_atom_species(ifa,ofs_running));
	EXPECT_DOUBLE_EQ(ucell->latvec.e11,4.27957);
	EXPECT_DOUBLE_EQ(ucell->latvec.e22,4.27957);
	EXPECT_DOUBLE_EQ(ucell->latvec.e33,4.27957);
	//mandatory preliminaries
	delete[] ucell->magnet.start_magnetization;
	ucell->magnet.start_magnetization = new double[ucell->ntype];
	ucell->read_atom_positions(ifa,ofs_running,ofs_warning);
	ofs_running.close();
	ofs_warning.close();
	ifa.close();
	remove("read_atom_positions.tmp");
	remove("read_atom_positions.warn");
}

TEST_F(UcellTest,ReadAtomPositionsCACXY)
{
	std::string fn = "./support/STRU_MgO_cacxy";
	std::ifstream ifa(fn.c_str());
	std::ofstream ofs_running;
	std::ofstream ofs_warning;
	ofs_running.open("read_atom_positions.tmp");
	ofs_warning.open("read_atom_positions.warn");
	//mandatory preliminaries
	ucell->ntype = 2;
	ucell->atoms = new Atom[ucell->ntype];
	ucell->set_atom_flag = true;
	GlobalV::test_pseudo_cell = 2;
	GlobalV::BASIS_TYPE = "lcao";
	GlobalV::deepks_setorb = true;
	GlobalV::NSPIN = 1;
	EXPECT_NO_THROW(ucell->read_atom_species(ifa,ofs_running));
	EXPECT_DOUBLE_EQ(ucell->latvec.e11,4.27957);
	EXPECT_DOUBLE_EQ(ucell->latvec.e22,4.27957);
	EXPECT_DOUBLE_EQ(ucell->latvec.e33,4.27957);
	//mandatory preliminaries
	delete[] ucell->magnet.start_magnetization;
	ucell->magnet.start_magnetization = new double[ucell->ntype];
	ucell->read_atom_positions(ifa,ofs_running,ofs_warning);
	ofs_running.close();
	ofs_warning.close();
	ifa.close();
	remove("read_atom_positions.tmp");
	remove("read_atom_positions.warn");
}

TEST_F(UcellTest,ReadAtomPositionsCACXZ)
{
	std::string fn = "./support/STRU_MgO_cacxz";
	std::ifstream ifa(fn.c_str());
	std::ofstream ofs_running;
	std::ofstream ofs_warning;
	ofs_running.open("read_atom_positions.tmp");
	ofs_warning.open("read_atom_positions.warn");
	//mandatory preliminaries
	ucell->ntype = 2;
	ucell->atoms = new Atom[ucell->ntype];
	ucell->set_atom_flag = true;
	GlobalV::test_pseudo_cell = 2;
	GlobalV::BASIS_TYPE = "lcao";
	GlobalV::deepks_setorb = true;
	GlobalV::NSPIN = 1;
	EXPECT_NO_THROW(ucell->read_atom_species(ifa,ofs_running));
	EXPECT_DOUBLE_EQ(ucell->latvec.e11,4.27957);
	EXPECT_DOUBLE_EQ(ucell->latvec.e22,4.27957);
	EXPECT_DOUBLE_EQ(ucell->latvec.e33,4.27957);
	//mandatory preliminaries
	delete[] ucell->magnet.start_magnetization;
	ucell->magnet.start_magnetization = new double[ucell->ntype];
	ucell->read_atom_positions(ifa,ofs_running,ofs_warning);
	ofs_running.close();
	ofs_warning.close();
	ifa.close();
	remove("read_atom_positions.tmp");
	remove("read_atom_positions.warn");
}

TEST_F(UcellTest,ReadAtomPositionsCACYZ)
{
	std::string fn = "./support/STRU_MgO_cacyz";
	std::ifstream ifa(fn.c_str());
	std::ofstream ofs_running;
	std::ofstream ofs_warning;
	ofs_running.open("read_atom_positions.tmp");
	ofs_warning.open("read_atom_positions.warn");
	//mandatory preliminaries
	ucell->ntype = 2;
	ucell->atoms = new Atom[ucell->ntype];
	ucell->set_atom_flag = true;
	GlobalV::test_pseudo_cell = 2;
	GlobalV::BASIS_TYPE = "lcao";
	GlobalV::deepks_setorb = true;
	GlobalV::NSPIN = 1;
	EXPECT_NO_THROW(ucell->read_atom_species(ifa,ofs_running));
	EXPECT_DOUBLE_EQ(ucell->latvec.e11,4.27957);
	EXPECT_DOUBLE_EQ(ucell->latvec.e22,4.27957);
	EXPECT_DOUBLE_EQ(ucell->latvec.e33,4.27957);
	//mandatory preliminaries
	delete[] ucell->magnet.start_magnetization;
	ucell->magnet.start_magnetization = new double[ucell->ntype];
	ucell->read_atom_positions(ifa,ofs_running,ofs_warning);
	ofs_running.close();
	ofs_warning.close();
	ifa.close();
	remove("read_atom_positions.tmp");
	remove("read_atom_positions.warn");
}

TEST_F(UcellTest,ReadAtomPositionsCACXYZ)
{
	std::string fn = "./support/STRU_MgO_cacxyz";
	std::ifstream ifa(fn.c_str());
	std::ofstream ofs_running;
	std::ofstream ofs_warning;
	ofs_running.open("read_atom_positions.tmp");
	ofs_warning.open("read_atom_positions.warn");
	//mandatory preliminaries
	ucell->ntype = 2;
	ucell->atoms = new Atom[ucell->ntype];
	ucell->set_atom_flag = true;
	GlobalV::test_pseudo_cell = 2;
	GlobalV::BASIS_TYPE = "lcao";
	GlobalV::deepks_setorb = true;
	GlobalV::NSPIN = 1;
	EXPECT_NO_THROW(ucell->read_atom_species(ifa,ofs_running));
	EXPECT_DOUBLE_EQ(ucell->latvec.e11,4.27957);
	EXPECT_DOUBLE_EQ(ucell->latvec.e22,4.27957);
	EXPECT_DOUBLE_EQ(ucell->latvec.e33,4.27957);
	//mandatory preliminaries
	delete[] ucell->magnet.start_magnetization;
	ucell->magnet.start_magnetization = new double[ucell->ntype];
	ucell->read_atom_positions(ifa,ofs_running,ofs_warning);
	ofs_running.close();
	ofs_warning.close();
	ifa.close();
	remove("read_atom_positions.tmp");
	remove("read_atom_positions.warn");
}

TEST_F(UcellTest,ReadAtomPositionsCAU)
{
	std::string fn = "./support/STRU_MgO_cau";
	std::ifstream ifa(fn.c_str());
	std::ofstream ofs_running;
	std::ofstream ofs_warning;
	ofs_running.open("read_atom_positions.tmp");
	ofs_warning.open("read_atom_positions.warn");
	//mandatory preliminaries
	ucell->ntype = 2;
	ucell->atoms = new Atom[ucell->ntype];
	ucell->set_atom_flag = true;
	GlobalV::test_pseudo_cell = 2;
	GlobalV::BASIS_TYPE = "lcao";
	GlobalV::deepks_setorb = true;
	GlobalV::NSPIN = 1;
	GlobalV::fixed_atoms = true;
	EXPECT_NO_THROW(ucell->read_atom_species(ifa,ofs_running));
	EXPECT_DOUBLE_EQ(ucell->latvec.e11,4.27957);
	EXPECT_DOUBLE_EQ(ucell->latvec.e22,4.27957);
	EXPECT_DOUBLE_EQ(ucell->latvec.e33,4.27957);
	//mandatory preliminaries
	delete[] ucell->magnet.start_magnetization;
	ucell->magnet.start_magnetization = new double[ucell->ntype];
	ucell->read_atom_positions(ifa,ofs_running,ofs_warning);
	ofs_running.close();
	ofs_warning.close();
	ifa.close();
	remove("read_atom_positions.tmp");
	remove("read_atom_positions.warn");
}

TEST_F(UcellTest,ReadAtomPositionsWarning1)
{
	std::string fn = "./support/STRU_MgO_WarningC1";
	std::ifstream ifa(fn.c_str());
	std::ofstream ofs_running;
	std::ofstream ofs_warning;
	ofs_running.open("read_atom_positions.tmp");
	ofs_warning.open("read_atom_positions.warn");
	//mandatory preliminaries
	ucell->ntype = 2;
	ucell->atoms = new Atom[ucell->ntype];
	ucell->set_atom_flag = true;
	GlobalV::test_pseudo_cell = 2;
	GlobalV::BASIS_TYPE = "lcao";
	GlobalV::deepks_setorb = true;
	EXPECT_NO_THROW(ucell->read_atom_species(ifa,ofs_running));
	EXPECT_DOUBLE_EQ(ucell->latvec.e11,4.27957);
	EXPECT_DOUBLE_EQ(ucell->latvec.e22,4.27957);
	EXPECT_DOUBLE_EQ(ucell->latvec.e33,4.27957);
	//mandatory preliminaries
	delete[] ucell->magnet.start_magnetization;
	ucell->magnet.start_magnetization = new double[ucell->ntype];
	EXPECT_NO_THROW(ucell->read_atom_positions(ifa,ofs_running,ofs_warning));
	ofs_running.close();
	ofs_warning.close();
	ifa.close();
	//check warning file
	std::ifstream ifs_tmp;
	ifs_tmp.open("read_atom_positions.warn");
	std::string str((std::istreambuf_iterator<char>(ifs_tmp)),std::istreambuf_iterator<char>());
	EXPECT_THAT(str, testing::HasSubstr("There are several options for you:"));
	EXPECT_THAT(str, testing::HasSubstr("Direct"));
	EXPECT_THAT(str, testing::HasSubstr("Cartesian_angstrom"));
	EXPECT_THAT(str, testing::HasSubstr("Cartesian_au"));
	EXPECT_THAT(str, testing::HasSubstr("Cartesian_angstrom_center_xy"));
	EXPECT_THAT(str, testing::HasSubstr("Cartesian_angstrom_center_xz"));
	EXPECT_THAT(str, testing::HasSubstr("Cartesian_angstrom_center_yz"));
	EXPECT_THAT(str, testing::HasSubstr("Cartesian_angstrom_center_xyz"));
	ifs_tmp.close();
	remove("read_atom_positions.tmp");
	remove("read_atom_positions.warn");
}

TEST_F(UcellTest,ReadAtomPositionsWarning2)
{
	std::string fn = "./support/STRU_MgO_WarningC2";
	std::ifstream ifa(fn.c_str());
	std::ofstream ofs_running;
	std::ofstream ofs_warning;
	ofs_running.open("read_atom_positions.tmp");
	ofs_warning.open("read_atom_positions.warn");
	//mandatory preliminaries
	ucell->ntype = 2;
	ucell->atoms = new Atom[ucell->ntype];
	ucell->set_atom_flag = true;
	GlobalV::test_pseudo_cell = 2;
	GlobalV::BASIS_TYPE = "lcao";
	GlobalV::deepks_setorb = true;
	EXPECT_NO_THROW(ucell->read_atom_species(ifa,ofs_running));
	EXPECT_DOUBLE_EQ(ucell->latvec.e11,4.27957);
	EXPECT_DOUBLE_EQ(ucell->latvec.e22,4.27957);
	EXPECT_DOUBLE_EQ(ucell->latvec.e33,4.27957);
	//mandatory preliminaries
	delete[] ucell->magnet.start_magnetization;
	ucell->magnet.start_magnetization = new double[ucell->ntype];
	EXPECT_NO_THROW(ucell->read_atom_positions(ifa,ofs_running,ofs_warning));
	ofs_running.close();
	ofs_warning.close();
	ifa.close();
	//check warning file
	std::ifstream ifs_tmp;
	ifs_tmp.open("read_atom_positions.warn");
	std::string str((std::istreambuf_iterator<char>(ifs_tmp)),std::istreambuf_iterator<char>());
	EXPECT_THAT(str, testing::HasSubstr("Label read from ATOMIC_POSITIONS is Mo"));
	EXPECT_THAT(str, testing::HasSubstr("Label from ATOMIC_SPECIES is Mg"));
	ifs_tmp.close();
	remove("read_atom_positions.tmp");
	remove("read_atom_positions.warn");
}

TEST_F(UcellTest,ReadAtomPositionsWarning3)
{
	std::string fn = "./support/STRU_MgO_WarningC3";
	std::ifstream ifa(fn.c_str());
	std::ofstream ofs_running;
	ofs_running.open("read_atom_positions.tmp");
	GlobalV::ofs_warning.open("read_atom_positions.warn");
	//mandatory preliminaries
	ucell->ntype = 2;
	ucell->atoms = new Atom[ucell->ntype];
	ucell->set_atom_flag = true;
	GlobalV::test_pseudo_cell = 2;
	GlobalV::BASIS_TYPE = "lcao";
	GlobalV::deepks_setorb = true;
	EXPECT_NO_THROW(ucell->read_atom_species(ifa,ofs_running));
	EXPECT_DOUBLE_EQ(ucell->latvec.e11,4.27957);
	EXPECT_DOUBLE_EQ(ucell->latvec.e22,4.27957);
	EXPECT_DOUBLE_EQ(ucell->latvec.e33,4.27957);
	//mandatory preliminaries
	delete[] ucell->magnet.start_magnetization;
	ucell->magnet.start_magnetization = new double[ucell->ntype];
	EXPECT_NO_THROW(ucell->read_atom_positions(ifa,ofs_running,GlobalV::ofs_warning));
	ofs_running.close();
	GlobalV::ofs_warning.close();
	ifa.close();
	//check warning file
	std::ifstream ifs_tmp;
	ifs_tmp.open("read_atom_positions.warn");
	std::string str((std::istreambuf_iterator<char>(ifs_tmp)),std::istreambuf_iterator<char>());
	EXPECT_THAT(str, testing::HasSubstr("read_atom_positions  warning :  atom number < 0."));
	ifs_tmp.close();
	remove("read_atom_positions.tmp");
	remove("read_atom_positions.warn");
}

TEST_F(UcellDeathTest,ReadAtomPositionsWarning4)
{
	std::string fn = "./support/STRU_MgO_WarningC4";
	std::ifstream ifa(fn.c_str());
	std::ofstream ofs_running;
	std::ofstream ofs_warning;
	ofs_running.open("read_atom_positions.tmp");
	ofs_warning.open("read_atom_positions.warn");
	//mandatory preliminaries
	ucell->ntype = 2;
	ucell->atoms = new Atom[ucell->ntype];
	ucell->set_atom_flag = true;
	GlobalV::test_pseudo_cell = 2;
	GlobalV::BASIS_TYPE = "lcao";
	GlobalV::deepks_setorb = true;
	EXPECT_NO_THROW(ucell->read_atom_species(ifa,ofs_running));
	EXPECT_DOUBLE_EQ(ucell->latvec.e11,4.27957);
	EXPECT_DOUBLE_EQ(ucell->latvec.e22,4.27957);
	EXPECT_DOUBLE_EQ(ucell->latvec.e33,4.27957);
	//mandatory preliminaries
	delete[] ucell->magnet.start_magnetization;
	ucell->magnet.start_magnetization = new double[ucell->ntype];
	testing::internal::CaptureStdout();
	EXPECT_EXIT(ucell->read_atom_positions(ifa,ofs_running,ofs_warning),
			::testing::ExitedWithCode(1),"");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("read_atom_positions, mismatch in atom number for atom type: Mg"));
	ofs_running.close();
	ofs_warning.close();
	ifa.close();
	remove("read_atom_positions.tmp");
	remove("read_atom_positions.warn");
}

TEST_F(UcellTest,ReadAtomPositionsWarning5)
{
	std::string fn = "./support/STRU_MgO";
	std::ifstream ifa(fn.c_str());
	std::ofstream ofs_running;
	ofs_running.open("read_atom_positions.tmp");
	GlobalV::ofs_warning.open("read_atom_positions.warn");
	//mandatory preliminaries
	ucell->ntype = 2;
	ucell->atoms = new Atom[ucell->ntype];
	ucell->set_atom_flag = true;
	GlobalV::test_pseudo_cell = 2;
	GlobalV::BASIS_TYPE = "lcao";
	GlobalV::deepks_setorb = true;
	GlobalV::CALCULATION = "md";
	GlobalV::ESOLVER_TYPE="arbitrary";
	EXPECT_NO_THROW(ucell->read_atom_species(ifa,ofs_running));
	EXPECT_DOUBLE_EQ(ucell->latvec.e11,4.27957);
	EXPECT_DOUBLE_EQ(ucell->latvec.e22,4.27957);
	EXPECT_DOUBLE_EQ(ucell->latvec.e33,4.27957);
	//mandatory preliminaries
	delete[] ucell->magnet.start_magnetization;
	ucell->magnet.start_magnetization = new double[ucell->ntype];
	EXPECT_NO_THROW(ucell->read_atom_positions(ifa,ofs_running,GlobalV::ofs_warning));
	ofs_running.close();
	GlobalV::ofs_warning.close();
	ifa.close();
	//check warning file
	std::ifstream ifs_tmp;
	ifs_tmp.open("read_atom_positions.warn");
	std::string str((std::istreambuf_iterator<char>(ifs_tmp)),std::istreambuf_iterator<char>());
	EXPECT_THAT(str, testing::HasSubstr("read_atoms  warning : no atom can move in MD!"));
	ifs_tmp.close();
	remove("read_atom_positions.tmp");
	remove("read_atom_positions.warn");
}
#endif
