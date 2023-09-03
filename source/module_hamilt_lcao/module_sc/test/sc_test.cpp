#include "../spin_constrain.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"

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

TEST_F(SpinConstrainTest, Nat)
{
    sc.set_nat(10);
    EXPECT_EQ(sc.get_nat(), 10);
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