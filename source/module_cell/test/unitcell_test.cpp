#include "gtest/gtest.h"
#include "module_base/mathzone.h"
#include "module_cell/unitcell.h"

InfoNonlocal::InfoNonlocal(){}
InfoNonlocal::~InfoNonlocal(){}
Magnetism::Magnetism(){}
Magnetism::~Magnetism(){}

/************************************************
 *  unit test of class UnitCell
 ***********************************************/

/**
 * - Tested Functions:
 *   - Constructor: UnitCell() and ~UnitCell()
 *
 */

//mock function
void LCAO_Orbitals::bcast_files(
	const int &ntype_in,
	const int &my_rank)
{
	return;
}

class UnitCellTest : public testing::Test
{
protected:
	std::unique_ptr<UnitCell> ucell{new UnitCell};
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
		}
		else if(fixed_axes_in[i] == "a")
		{
			EXPECT_EQ(ucell->lc[0],0);
			EXPECT_EQ(ucell->lc[1],1);
			EXPECT_EQ(ucell->lc[2],1);
		}
		else if(fixed_axes_in[i] == "b")
		{
			EXPECT_EQ(ucell->lc[0],1);
			EXPECT_EQ(ucell->lc[1],0);
			EXPECT_EQ(ucell->lc[2],1);
		}
		else if(fixed_axes_in[i] == "c")
		{
			EXPECT_EQ(ucell->lc[0],1);
			EXPECT_EQ(ucell->lc[1],1);
			EXPECT_EQ(ucell->lc[2],0);
		}
		else if(fixed_axes_in[i] == "ab")
		{
			EXPECT_EQ(ucell->lc[0],0);
			EXPECT_EQ(ucell->lc[1],1);
			EXPECT_EQ(ucell->lc[2],1);
		}
		else if(fixed_axes_in[i] == "ac")
		{
			EXPECT_EQ(ucell->lc[0],0);
			EXPECT_EQ(ucell->lc[1],1);
			EXPECT_EQ(ucell->lc[2],0);
		}
		else if(fixed_axes_in[i] == "bc")
		{
			EXPECT_EQ(ucell->lc[0],1);
			EXPECT_EQ(ucell->lc[1],0);
			EXPECT_EQ(ucell->lc[2],0);
		}
		else if(fixed_axes_in[i] == "abc")
		{
			EXPECT_EQ(ucell->lc[0],0);
			EXPECT_EQ(ucell->lc[1],0);
			EXPECT_EQ(ucell->lc[2],0);
		}
	}
}

TEST_F(UnitCellTest,SetAtom)
{
	ucell->ntype = 1;
	ucell->atoms = new Atom[ucell->ntype];
	ucell->set_atom_flag = true;
	delete[] ucell->atom_label;
	delete[] ucell->atom_mass;
	delete[] ucell->pseudo_fn;
	delete[] ucell->pseudo_type;
	delete[] ucell->orbital_fn;
	ucell->atom_label = new std::string[ucell->ntype];
	ucell->atom_mass = new double[ucell->ntype];
	ucell->pseudo_fn = new std::string[ucell->ntype];
	ucell->pseudo_type = new std::string[ucell->ntype];
	ucell->orbital_fn = new std::string[ucell->ntype];
	ucell->atom_label[0] = "C";
	ucell->atom_mass[0] = 12.0;
	ucell->pseudo_fn[0] = "C.upf";
	ucell->pseudo_type[0] = "upf201";
	ucell->orbital_fn[0] = "C.orb";
	ucell->lat0 = 1.8897261254578281;
	ucell->tpiba = ModuleBase::TWO_PI/ucell->lat0;
	ucell->tpiba2 = ucell->tpiba * ucell->tpiba;
	ucell->latvec.e11 = 10.0; ucell->latvec.e12 = 0.0; ucell->latvec.e13 = 0.0;
	ucell->latvec.e21 = 0.0; ucell->latvec.e22 = 10.0; ucell->latvec.e23 = 0.0;
	ucell->latvec.e31 = 0.0; ucell->latvec.e32 = 0.0; ucell->latvec.e33 = 10.0;
	ucell->GT = ucell->latvec.Inverse();
	ucell->G = ucell->GT.Transpose();
	ucell->a1.x = ucell->latvec.e11; ucell->a1.y = ucell->latvec.e12; ucell->a1.z = ucell->latvec.e13;
	ucell->a2.x = ucell->latvec.e21; ucell->a2.y = ucell->latvec.e22; ucell->a2.z = ucell->latvec.e23;
	ucell->a3.x = ucell->latvec.e31; ucell->a3.y = ucell->latvec.e32; ucell->a3.z = ucell->latvec.e33;
	ucell->Coordinate = "Direct";
	ucell->atoms[0].label = "C";
	//nw
	ucell->atoms[0].nw = 0;
	ucell->atoms[0].nwl = 2;
	ucell->lmaxmax = 2;
	delete[] ucell->atoms[0].l_nchi;
	ucell->atoms[0].l_nchi = new int[ ucell->atoms[0].nwl+1];
	for(int L=0; L<ucell->atoms[0].nwl+1; L++)
	{
		ucell->atoms[0].l_nchi[L] = 1;
		ucell->atoms[0].nw += (2*L + 1) * ucell->atoms[0].l_nchi[L];
	}
}
