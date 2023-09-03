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
};

TEST_F(SpinConstrainTest, Nat)
{
    SpinConstrain& sc = SpinConstrain::getInstance();
    sc.set_nat(10);
    EXPECT_EQ(sc.get_nat(), 10);
}

/*
TEST_F(SpinConstrainTest, PaserScJsonFile)
{
	std::map<std::string, std::vector<Input_Conv::ScElementData>> data1;
    parseScJsonFile("./support/sc.json", data1);
	EXPECT_EQ(data1.size(), 2);
	for (const auto& sc : data1)
	{
        const std::string& element_name = sc.first;
        const std::vector<Input_Conv::ScElementData>& sc_atoms = sc.second;
		if (element_name == "Si")
		{
			EXPECT_EQ(sc_atoms.size(), 2);
		}
		else if (element_name == "O")
		{
			EXPECT_EQ(sc_atoms.size(), 1);
		}
        for (const Input_Conv::ScElementData& sc_data : sc_atoms) {
			if (element_name == "Si" & sc_data.index == 0)
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
    std::map<std::string, std::vector<Input_Conv::ScElementData>> data2;
    parseScJsonFile("./support/sc_v2.json", data2);
	EXPECT_EQ(data2.size(), 1);
	for (const auto& sc : data2)
	{
        const std::string& element_name = sc.first;
        const std::vector<Input_Conv::ScElementData>& sc_atoms = sc.second;
		EXPECT_EQ(element_name, "Si");
		EXPECT_EQ(sc_atoms.size(), 2);
        for (const Input_Conv::ScElementData& sc_data : sc_atoms) {
			if (element_name == "Si" & sc_data.index == 0)
			{
				EXPECT_DOUBLE_EQ(sc_data.lambda[0],0.1);
				EXPECT_DOUBLE_EQ(sc_data.lambda[1],0.1);
				EXPECT_DOUBLE_EQ(sc_data.lambda[2],0.2);
				EXPECT_DOUBLE_EQ(sc_data.sc_spin_val,1.0);
				EXPECT_DOUBLE_EQ(sc_data.sc_spin_angle1,30.0);
				EXPECT_DOUBLE_EQ(sc_data.sc_spin_angle2,45.0);
			}
        }
	}
}
*/