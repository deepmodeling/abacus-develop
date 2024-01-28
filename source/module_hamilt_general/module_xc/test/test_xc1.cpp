#include "../xc_functional.h"
#include "gtest/gtest.h"
#include "../exx_info.h"

/************************************************
*  unit test of set_xc_type
***********************************************/

// For more information of the functions, check the comment of xc_functional.h
// the functionals are not tested because they all use libxc
// so only set_xc_type is called

namespace ModuleBase
{
    void WARNING_QUIT(const std::string &file,const std::string &description) {return ;}
}

namespace GlobalV
{
    std::string BASIS_TYPE = "";
    bool CAL_STRESS = 0;
    int CAL_FORCE = 0;
    int NSPIN = 1;
    double XC_TEMPERATURE;
}

namespace GlobalC
{
	Exx_Info exx_info;
}

class XCTest_HSE : public testing::Test
{
    protected:
        std::vector<double> e_lda, v_lda;
        std::vector<double> e_gga, v1_gga, v2_gga;

        void SetUp()
        {
            XC_Functional::set_xc_type("HSE");
            XC_Functional::get_hybrid_alpha(0.5);
        }
};

TEST_F(XCTest_HSE, set_xc_type)
{
    EXPECT_EQ(XC_Functional::get_func_type(),4);
}

class XCTest_SCAN0 : public testing::Test
{
    protected:
        void SetUp()
        {
            XC_Functional::set_xc_type("SCAN0");
            XC_Functional::get_hybrid_alpha(0.5);
        }
};

TEST_F(XCTest_SCAN0, set_xc_type)
{
    EXPECT_EQ(XC_Functional::get_func_type(),5);
}

class XCTest_KSDT : public testing::Test
{
    protected:
        void SetUp()
        {
            XC_Functional::set_xc_type("XC_LDA_XC_KSDT");
            GlobalV::XC_TEMPERATURE = 0.5;
        }
};

TEST_F(XCTest_KSDT, set_xc_type)
{
    EXPECT_EQ(XC_Functional::get_func_type(),1);
}

class XCTest_KT2 : public testing::Test
{
    protected:
        void SetUp()
        {
            XC_Functional::set_xc_type("GGA_XC_KT2");
        }
};

TEST_F(XCTest_KT2, set_xc_type)
{
    EXPECT_EQ(XC_Functional::get_func_type(),2);
}

class XCTest_R2SCAN : public testing::Test
{
    protected:
        void SetUp()
        {
            XC_Functional::set_xc_type("MGGA_X_R2SCAN+MGGA_C_R2SCAN");
        }
};

TEST_F(XCTest_R2SCAN, set_xc_type)
{
    EXPECT_EQ(XC_Functional::get_func_type(),3);
}

class XCTest_LB07 : public testing::Test
{
    protected:
        void SetUp()
        {
            XC_Functional::set_xc_type("HYB_GGA_XC_LB07");
        }
};

TEST_F(XCTest_LB07, set_xc_type)
{
    EXPECT_EQ(XC_Functional::get_func_type(),4);
}

class XCTest_BMK : public testing::Test
{
    protected:
        void SetUp()
        {
            XC_Functional::set_xc_type("HYB_MGGA_X_BMK");
        }
};

TEST_F(XCTest_BMK, set_xc_type)
{
    EXPECT_EQ(XC_Functional::get_func_type(),5);
}

class XCTest_HF : public testing::Test
{
    protected:
        void SetUp()
        {
            XC_Functional::set_xc_type("HF");
        }
};

TEST_F(XCTest_HF, set_xc_type)
{
    EXPECT_EQ(XC_Functional::get_func_type(),4);
}