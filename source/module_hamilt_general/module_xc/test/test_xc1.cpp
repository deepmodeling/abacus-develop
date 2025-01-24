#include "gtest/gtest.h"
#include "xctest.h"
#include "../xc_functional.h"
#include "../exx_info.h"

/************************************************
*  unit test of set_xc_type
***********************************************/

// For more information of the functions, check the comment of xc_functional.h
// the functionals are not tested because they all use libxc
// so only set_xc_type is called

namespace ModuleBase
{
    void WARNING_QUIT(const std::string &file,const std::string &description) {exit(1);}
}

namespace GlobalV
{
    std::string BASIS_TYPE = "";
    bool CAL_STRESS = 0;
    int CAL_FORCE = 0;
    int NSPIN = 1;
}

namespace GlobalC
{
	Exx_Info exx_info;
}

class XCTest_HSE : public XCTest
{
    protected:
        std::vector<double> e_lda, v_lda;
        std::vector<double> e_gga, v1_gga, v2_gga;

        void SetUp()
        {
            XC_Functional::set_xc_type("HSE");
            XC_Functional::set_hybrid_alpha(0.5);
        }
};

TEST_F(XCTest_HSE, set_xc_type)
{
    EXPECT_EQ(XC_Functional::get_func_type(),4);
}

class XCTest_SCAN0 : public XCTest
{
    protected:
        void SetUp()
        {
            XC_Functional::set_xc_type("SCAN0");
            XC_Functional::set_hybrid_alpha(0.5);
        }
};

TEST_F(XCTest_SCAN0, set_xc_type)
{
    EXPECT_EQ(XC_Functional::get_func_type(),5);
    EXPECT_TRUE(XC_Functional::get_ked_flag());
}

class XCTest_KSDT : public XCTest
{
    protected:
        void SetUp()
        {
            XC_Functional::set_xc_type("XC_LDA_XC_KSDT");
        }
};

TEST_F(XCTest_KSDT, set_xc_type)
{
    EXPECT_EQ(XC_Functional::get_func_type(),1);
}

class XCTest_KT2 : public XCTest
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

class XCTest_R2SCAN : public XCTest
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
    EXPECT_TRUE(XC_Functional::get_ked_flag());
}

class XCTest_LB07 : public XCTest
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

class XCTest_BMK : public XCTest
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
    EXPECT_TRUE(XC_Functional::get_ked_flag());
}

class XCTest_HF : public XCTest
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