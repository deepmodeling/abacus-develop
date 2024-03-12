#include "../spin_constrain.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <algorithm>

/************************************************
 *  unit test of functions in class SpinConstrain
 ***********************************************/

/**
 * - Tested functions:
 *  - SpinConstrain::Set_ScData_From_Json()
 *     set the map from element index to ScAtomData from json file
 *  - SpinConstrain::get_ScData()
 *     get the map from element index to ScAtomData
 */

K_Vectors::K_Vectors(){}
K_Vectors::~K_Vectors(){}

template <typename T>
class SpinConstrainTest : public testing::Test
{
  protected:
    SpinConstrain<T, psi::DEVICE_CPU>& sc = SpinConstrain<T, psi::DEVICE_CPU>::getScInstance();
};

using MyTypes = ::testing::Types<double, std::complex<double>>;
TYPED_TEST_SUITE(SpinConstrainTest, MyTypes);

TYPED_TEST(SpinConstrainTest, ScDataFormat1)
{
    this->sc.Set_ScData_From_Json("./support/sc_f1.json");
    EXPECT_EQ(this->sc.get_ScData().size(), 2);
    for (const auto& sc_elem: this->sc.get_ScData())
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

TYPED_TEST(SpinConstrainTest, ScDataFormat2)
{
    this->sc.Set_ScData_From_Json("./support/sc_f2.json");
    EXPECT_EQ(this->sc.get_ScData().size(), 1);
    for (const auto& sc_elem: this->sc.get_ScData())
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
    EXPECT_DOUBLE_EQ(this->sc.get_decay_grad(1), 0.9);
}

TYPED_TEST(SpinConstrainTest, ScDataFormat3)
{
    this->sc.Set_ScData_From_Json("./support/sc_f3.json");
    EXPECT_EQ(this->sc.get_ScData().size(), 1);
    for (const auto& sc_elem: this->sc.get_ScData())
    {
        const int& it = sc_elem.first;
        const std::vector<ScAtomData>& sc_atoms = sc_elem.second;
		EXPECT_EQ(sc_atoms.size(), 2);
		EXPECT_EQ(it, 1);
        for (const ScAtomData& sc_data : sc_atoms) {
			if (it == 1 & sc_data.index == 4)
			{
				EXPECT_DOUBLE_EQ(sc_data.lambda[0],0.0);
				EXPECT_DOUBLE_EQ(sc_data.lambda[1],0.0);
				EXPECT_DOUBLE_EQ(sc_data.lambda[2],0.2);
				EXPECT_DOUBLE_EQ(sc_data.target_mag[0],0.0);
                EXPECT_DOUBLE_EQ(sc_data.target_mag[1],0.0);
                EXPECT_DOUBLE_EQ(sc_data.target_mag[2],10.0);
                EXPECT_EQ(sc_data.mag_type, 0);
            }
        }
	}
    EXPECT_DOUBLE_EQ(this->sc.get_decay_grad(1), 0.0);
}

TYPED_TEST(SpinConstrainTest, ScDataWarning)
{
    testing::internal::CaptureStdout();
    EXPECT_EXIT(this->sc.Set_ScData_From_Json("./support/sc_f4.json"), ::testing::ExitedWithCode(0), "");
    std::string output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("Error opening sc_file"));
}