#include "../spin_constrain.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <algorithm>

/************************************************
 *  unit test of functions in class SpinConstrain
 ***********************************************/

/**
 * - Tested functions:
 *  - SpinConstrain::getScInstance()
 *      get the instance of SpinConstrain
 *  - SpinConstrain::clear_atomCounts()
 *      which is a map from element index to atom number
 *  - SpinConstrain::set_atomCounts()
 *     set the map from element index to atom number
 *  - SpinConstrain::get_atomCounts()
 *     get the map from element index to atom number
 *  - SpinConstrain::get_nat()
 *     get the total number of atoms
 *  - SpinConstrain::get_iat()
 *     get the atom index from (itype, atom_index)
 *  - SpinConstrain::clear_orbitalCounts()
 *     clear the map from element index to orbital number
 *  - SpinConstrain::set_orbitalCounts()
 *     set the map from element index to orbital number
 *  - SpinConstrain::get_orbitalCounts()
 *     get the map from element index to orbital number
 *  - SpinConstrain::get_nw()
 *     get the total number of orbitals
 *  - SpinConstrain::set_npol()
 *     set the number of npol, which is the number of spin components
 *  - SpinConstrain::get_npol()
 *     get the number of npol, which is the number of spin components
 *  - SpinConstrain::get_iwt()
 *     get the index of orbital with spin component from (itype, iat, orbital_index)
 */

K_Vectors::K_Vectors(){}
K_Vectors::~K_Vectors(){}

class SpinConstrainTest : public testing::Test
{
  protected:
  	SpinConstrain<std::complex<double>, psi::DEVICE_CPU>& sc = SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::getScInstance();
};

TEST_F(SpinConstrainTest, AtomCounts)
{
	std::map<int, int> atomCounts = {{0,5},{1,10}};
	sc.clear_atomCounts();
	sc.set_atomCounts(atomCounts);
	std::map<int, int> atomCounts2 = sc.get_atomCounts();
	int ntype = atomCounts2.size();
	EXPECT_EQ(ntype, 2);
    int nat = sc.get_nat();
	EXPECT_EQ(nat, 15);
    EXPECT_EQ(sc.get_iat(1, 4), 9); // atom_index starts from 0
}

TEST_F(SpinConstrainTest, OrbitalCounts)
{
	std::map<int, int> orbitalCounts = {{0,5},{1,10}};
	std::map<int, int> atomCounts = {{0,1},{1,2}};
	sc.clear_atomCounts();
	sc.clear_orbitalCounts();
	sc.set_atomCounts(atomCounts);
	sc.set_orbitalCounts(orbitalCounts);
	std::map<int, int> orbitalCounts2 = sc.get_orbitalCounts();
	int ntype = orbitalCounts2.size();
	EXPECT_EQ(ntype, 2);
	EXPECT_EQ(sc.get_nw(), 25);
	sc.set_npol(2);
	EXPECT_EQ(sc.get_npol(),2);
	EXPECT_EQ(sc.get_nw(), 50); // npol = 2
	sc.set_npol(1);
	EXPECT_EQ(sc.get_npol(),1);
	EXPECT_EQ(sc.get_iwt(1,1,2), 17);
	sc.set_npol(2);
	EXPECT_EQ(sc.get_iwt(1,1,2), 32); // npol = 2
}

TEST_F(SpinConstrainTest, SetScLambdaMagConstrain)
{
	sc.clear_ScData();
	sc.Set_ScData_From_Json("./support/sc_f1.json");
	std::map<int, int> atomCounts = {{0,5},{1,10}};
	sc.clear_atomCounts();
	sc.set_atomCounts(atomCounts);
	int nat = sc.get_nat();
	sc.set_sc_lambda();
	sc.set_target_mag();
    sc.set_constrain();
    std::vector<ModuleBase::Vector3<double>> sc_lambda = sc.get_sc_lambda();
    std::vector<ModuleBase::Vector3<double>> target_mag = sc.get_target_mag();
	std::vector<ModuleBase::Vector3<int>> constrain = sc.get_constrain();
    sc.set_sc_lambda(sc_lambda.data(), nat);
	sc.set_target_mag(target_mag.data(), nat);
	sc.set_constrain(constrain.data(), nat);
    EXPECT_EQ(sc_lambda.size(), sc.get_nat());
    for (const auto& sc_elem : sc.get_ScData())
	{
		int itype = sc_elem.first;
        const std::vector<ScAtomData>& sc_atoms = sc_elem.second;
        for (const ScAtomData& sc_data : sc_atoms) {
			int index = sc_data.index;
			int iat = sc.get_iat(itype, index);
			EXPECT_DOUBLE_EQ(sc_data.lambda[0]*sc.meV_to_Ry,sc_lambda[iat].x);
			EXPECT_DOUBLE_EQ(sc_data.lambda[1]*sc.meV_to_Ry,sc_lambda[iat].y);
			EXPECT_DOUBLE_EQ(sc_data.lambda[2]*sc.meV_to_Ry,sc_lambda[iat].z);
			EXPECT_DOUBLE_EQ(sc_data.target_mag[0],target_mag[iat].x);
			EXPECT_DOUBLE_EQ(sc_data.target_mag[1],target_mag[iat].y);
			EXPECT_DOUBLE_EQ(sc_data.target_mag[2],target_mag[iat].z);
			EXPECT_EQ(sc_data.constrain[0],constrain[iat].x);
			EXPECT_EQ(sc_data.constrain[1],constrain[iat].y);
			EXPECT_EQ(sc_data.constrain[2],constrain[iat].z);
		}
	}
	for (int iat = 0; iat < sc.get_nat(); iat++)
	{
		if (! (iat == 1 || iat ==5 || iat ==9))
		{
			EXPECT_DOUBLE_EQ(sc_lambda[iat].x,0.0);
			EXPECT_DOUBLE_EQ(sc_lambda[iat].y,0.0);
			EXPECT_DOUBLE_EQ(sc_lambda[iat].z,0.0);
			EXPECT_DOUBLE_EQ(target_mag[iat].x,0.0);
			EXPECT_DOUBLE_EQ(target_mag[iat].y,0.0);
			EXPECT_DOUBLE_EQ(target_mag[iat].z,0.0);
		}
	}
}

TEST_F(SpinConstrainTest, CalEscon)
{
    sc.zero_Mi();
	sc.clear_ScData();
	sc.Set_ScData_From_Json("./support/sc_f1.json");
	std::map<int, int> atomCounts = {{0,5},{1,10}};
	sc.clear_atomCounts();
	sc.set_atomCounts(atomCounts);
	int nat = sc.get_nat();
	sc.set_sc_lambda();
	double escon = sc.cal_escon();
    double escon1 = sc.get_escon();
	EXPECT_DOUBLE_EQ(escon, escon1);
    EXPECT_DOUBLE_EQ(escon1, 0.0);
}

TEST_F(SpinConstrainTest, NSPIN)
{
    sc.set_nspin(4);
	int nspin = sc.get_nspin();
	EXPECT_EQ(nspin, 4);
}