#include "../spin_constrain.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <algorithm>

/************************************************
 *  unit test of functions in class SpinConstrain
 ***********************************************/

/**
 * 
 * 
 * 
 */

class SpinConstrainTest : public testing::Test
{
  protected:
	SpinConstrain& sc = SpinConstrain::getInstance();
};

TEST_F(SpinConstrainTest, AtomCounts)
{
	std::map<int, int> atomCounts = {{0,5},{1,10}};
	sc.set_atomCounts(atomCounts);
	std::map<int, int> atomCounts2 = sc.get_atomCounts();
	int ntype = atomCounts2.size();
	EXPECT_EQ(ntype, 2);
    int nat = sc.get_nat();
	EXPECT_EQ(nat, 15);

}

TEST_F(SpinConstrainTest, ScDataFormat1)
{
	sc.clear_ScData();
	sc.Set_ScData_From_Json("./support/sc_f1.json");
	EXPECT_EQ(sc.get_ScData().size(), 2);
	for (const auto& sc_elem : sc.get_ScData())
	{
        const int& it = sc_elem.first;
        const std::vector<ScAtomData>& sc_atoms = sc_elem.second;
		if (it == 0)
		{
			EXPECT_EQ(sc_atoms.size(), 1);
		}
		else if (it == 1)
		{
			EXPECT_EQ(sc_atoms.size(), 2);
		}
        for (const ScAtomData& sc_data : sc_atoms) {
			if (it == 1 & sc_data.index == 0)
			{
				EXPECT_DOUBLE_EQ(sc_data.lambda[0],0.1);
				EXPECT_DOUBLE_EQ(sc_data.lambda[1],0.1);
				EXPECT_DOUBLE_EQ(sc_data.lambda[2],0.2);
				EXPECT_DOUBLE_EQ(sc_data.sc_mag[0],1.0);
				EXPECT_DOUBLE_EQ(sc_data.sc_mag[1],2.0);
				EXPECT_DOUBLE_EQ(sc_data.sc_mag[2],3.0);
			}
        }
	}
}

TEST_F(SpinConstrainTest, ScDataFormat2)
{
	sc.clear_ScData();
	sc.Set_ScData_From_Json("./support/sc_f2.json");
	EXPECT_EQ(sc.get_ScData().size(), 1);
	for (const auto& sc_elem : sc.get_ScData())
	{
        const int& it = sc_elem.first;
        const std::vector<ScAtomData>& sc_atoms = sc_elem.second;
		EXPECT_EQ(sc_atoms.size(), 2);
		EXPECT_EQ(it, 1);
        for (const ScAtomData& sc_data : sc_atoms) {
			if (it == 1 & sc_data.index == 4)
			{
				EXPECT_DOUBLE_EQ(sc_data.lambda[0],0.2);
				EXPECT_DOUBLE_EQ(sc_data.lambda[1],0.4);
				EXPECT_DOUBLE_EQ(sc_data.lambda[2],0.5);
				EXPECT_DOUBLE_EQ(sc_data.sc_spin_val,1.5);
				EXPECT_DOUBLE_EQ(sc_data.sc_spin_angle1,60.0);
				EXPECT_DOUBLE_EQ(sc_data.sc_spin_angle2,90.0);
			}
        }
	}
}

TEST_F(SpinConstrainTest, GetScLambdaAndMag)
{
	sc.clear_ScData();
	sc.Set_ScData_From_Json("./support/sc_f1.json");
	std::map<int, int> atomCounts = {{0,5},{1,10}};
	sc.clear_atomCounts();
	sc.set_atomCounts(atomCounts);
	std::vector<ModuleBase::Vector3<double>> sc_lambda = sc.get_sc_lambda();
	std::vector<ModuleBase::Vector3<double>> sc_mag = sc.get_sc_mag();
	EXPECT_EQ(sc_lambda.size(), sc.get_nat());
	for (const auto& sc_elem : sc.get_ScData())
	{
		int itype = sc_elem.first;
        const std::vector<ScAtomData>& sc_atoms = sc_elem.second;
        for (const ScAtomData& sc_data : sc_atoms) {
			int index = sc_data.index;
			int iat = sc.get_iat(itype, index);
			EXPECT_DOUBLE_EQ(sc_data.lambda[0],sc_lambda[iat].x);
			EXPECT_DOUBLE_EQ(sc_data.lambda[1],sc_lambda[iat].y);
			EXPECT_DOUBLE_EQ(sc_data.lambda[2],sc_lambda[iat].z);
			EXPECT_DOUBLE_EQ(sc_data.sc_mag[0],sc_mag[iat].x);
			EXPECT_DOUBLE_EQ(sc_data.sc_mag[1],sc_mag[iat].y);
			EXPECT_DOUBLE_EQ(sc_data.sc_mag[2],sc_mag[iat].z);
		}
	}
	for (int iat = 0; iat < sc.get_nat(); iat++)
	{
		if (! (iat == 1 || iat ==5 || iat ==9))
		{
			EXPECT_DOUBLE_EQ(sc_lambda[iat].x,0.0);
			EXPECT_DOUBLE_EQ(sc_lambda[iat].y,0.0);
			EXPECT_DOUBLE_EQ(sc_lambda[iat].z,0.0);
			EXPECT_DOUBLE_EQ(sc_mag[iat].x,0.0);
			EXPECT_DOUBLE_EQ(sc_mag[iat].y,0.0);
			EXPECT_DOUBLE_EQ(sc_mag[iat].z,0.0);
		}
	}
}