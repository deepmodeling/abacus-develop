#include "../spin_constrain.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <algorithm>

/************************************************
 *  unit test of functions in class SpinConstrain
 ***********************************************/

/**
 * - Tested functions:
 *  - SpinConstrain::clear_ScData()
 *     clear the map from element index to ScAtomData
 *  - SpinConstrain::Set_ScData_From_Json()
 *     set the map from element index to ScAtomData from json file
 *  - SpinConstrain::get_ScData()
 *     get the map from element index to ScAtomData
 */

K_Vectors::K_Vectors(){}
K_Vectors::~K_Vectors(){}

class SpinConstrainTest : public testing::Test
{
  protected:
  	SpinConstrain<std::complex<double>, psi::DEVICE_CPU>& sc = SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::getScInstance();
};

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
				EXPECT_DOUBLE_EQ(sc_data.target_mag[0],1.0);
				EXPECT_DOUBLE_EQ(sc_data.target_mag[1],2.0);
				EXPECT_DOUBLE_EQ(sc_data.target_mag[2],3.0);
				EXPECT_EQ(sc_data.constrain[0], 0);
				EXPECT_EQ(sc_data.constrain[1], 0);
				EXPECT_EQ(sc_data.constrain[2], 1);
                EXPECT_EQ(sc_data.mag_type, 0);
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
				EXPECT_DOUBLE_EQ(sc_data.target_mag_val,1.5);
                EXPECT_DOUBLE_EQ(sc_data.target_mag_angle1, 90.0);
                EXPECT_DOUBLE_EQ(sc_data.target_mag_angle2,90.0);
                EXPECT_EQ(sc_data.mag_type, 1);
            }
        }
	}
    EXPECT_DOUBLE_EQ(sc.get_decay_grad(1),0.9);
}

TEST_F(SpinConstrainTest, ScDataWarning)
{
	sc.clear_ScData();
	testing::internal::CaptureStdout();
	EXPECT_EXIT(sc.Set_ScData_From_Json("./support/sc_f3.json"), ::testing::ExitedWithCode(0), "");
	std::string output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("Error opening sc_file"));
}