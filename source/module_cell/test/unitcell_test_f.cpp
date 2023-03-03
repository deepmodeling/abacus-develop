#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "module_base/mathzone.h"
#include "module_cell/unitcell.h"

#ifdef __LCAO
InfoNonlocal::InfoNonlocal(){}
InfoNonlocal::~InfoNonlocal(){}
#endif
Magnetism::Magnetism(){}
Magnetism::~Magnetism(){}

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

class UnitCellTest : public testing::Test
{
protected:
	std::unique_ptr<UnitCell> ucell{new UnitCell};
	std::string output;
};

TEST_F(UnitCellTest,Constructor)
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

TEST_F(UnitCellTest,Setup)
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

TEST_F(UnitCellTest,SetupWarningQuit1)
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

TEST_F(UnitCellTest,SetupWarningQuit2)
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

TEST_F(UnitCellTest,RemakeCell)
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

TEST_F(UnitCellTest,RemakeCellWarnings)
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


TEST_F(UnitCellTest,JudgeParallel)
{
	ModuleBase::Vector3<double> b(1.0,1.0,1.0);
	double a[3] = {1.0,1.0,1.0};
	EXPECT_TRUE(ucell->judge_parallel(a,b));
}
