#include "../xc_functional.h"
#include "gtest/gtest.h"
#include "../exx_info.h"

/************************************************
*  unit test of functionals
***********************************************/

// For more information of the functions, check the comment of xc_functional.h
// mGGA Functional (SCAN) is tested under nspin = 1

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
    double XC_TEMPERATURE;
}

namespace GlobalC
{
	Exx_Info exx_info;
}

class XCTest_SCAN : public testing::Test
{
    protected:
        std::vector<double> e_,v1_,v2_,v3_;

        void SetUp()
        {
            XC_Functional::set_xc_type("SCAN");
            std::vector<double> rho  = {0.17E+01, 0.17E+01, 0.15E+01, 0.88E-01, 0.18E+04};
            std::vector<double> grho = {0.81E-11, 0.17E+01, 0.36E+02, 0.87E-01, 0.55E+00};
            std::vector<double> tau  = {0,0.02403590412,0.01672229351,0.01340429824,0.01141731056};

            for(int i=0;i<5;i++)
            {
                double e,v,v1,v2,v3;
                XC_Functional::tau_xc(rho[i],grho[i],tau[i],e,v1,v2,v3);
                e_.push_back(e);
                v1_.push_back(v1);
                v2_.push_back(v2);
                v3_.push_back(v3);
            }
        }
};

TEST_F(XCTest_SCAN, set_xc_type)
{
    EXPECT_EQ(XC_Functional::get_func_type(),3);
    std::vector<double> e_ref  = {-1.802228724,-1.802210367,-1.526620682,-0.03583190143,-19035.36827};
    std::vector<double> v1_ref = {-1.40580492,-1.405589907,-1.349177581,-0.5343722373,-14.09114168};
    std::vector<double> v2_ref = {-0.002302953645,-0.002306239202,-0.002688619206,-0.07677599434,-3.0633237e-07};
    std::vector<double> v3_ref = {0.01642445001,0.01644552589,0.01695685344,0.03307669533,0.002214113134};

    for (int i = 0;i<4;++i)
    {
        EXPECT_NEAR(e_[i],e_ref[i],1.0e-8);
        EXPECT_NEAR(v1_[i],v1_ref[i],1.0e-8);
        EXPECT_NEAR(v2_[i],v2_ref[i],1.0e-8);
        EXPECT_NEAR(v3_[i],v3_ref[i],1.0e-8);
    }
}