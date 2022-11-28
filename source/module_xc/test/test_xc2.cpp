#include "../xc_functional.h"
#include "gtest/gtest.h"
#include "../exx_info.h"

namespace ModuleBase
{
    void WARNING_QUIT(const std::string &file,const std::string &description) {return ;}
}

namespace GlobalV
{
    std::string BASIS_TYPE = "";
    bool CAL_STRESS = 0;
    int CAL_FORCE = 0;
    int NSPIN = 2;
    double XC_TEMPERATURE;
}

namespace GlobalC
{
	Exx_Info exx_info;
}

class XCTest_SLATER1_SPN : public testing::Test
{
    protected:
        std::vector<double> e_lda, v1_lda, v2_lda;

        void SetUp()
        {
            std::vector<double> rho  = {0.17E+01, 0.17E+01, 0.15E+01, 0.88E-01, 0.18E+04};
            std::vector<double> zeta = {0.0, 0.2, 0.5, 0.8, 1.0};

            for(int i=0;i<5;i++)
            {
                double e,v1,v2;
                XC_Functional::slater1_spin(rho[i],zeta[i],e,v1,v2);
                e_lda.push_back(e);
                v1_lda.push_back(v1);
                v2_lda.push_back(v2);
            }          
        }
};

TEST_F(XCTest_SLATER1_SPN, set_xc_type)
{
    std::vector<double> e_lda_ref  = {-1.32218621089,-1.33398308443,-1.34039342564,-0.568290937412,-16.9789364717};
    std::vector<double> v1_lda_ref  = {-1.76291494786,-1.87337667608,-1.93557153111,-0.799220801445,-22.6385819622};
    std::vector<double> v2_lda_ref  = {-1.76291494786,-1.63654526731,-1.34205034341,-0.384225285821,0};

    for (int i = 0;i<5;++i)
    {
        EXPECT_NEAR(e_lda[i],e_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v1_lda[i],v1_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v2_lda[i],v2_lda_ref[i],1.0e-8);
    }
}

class XCTest_SLATER_RXC_SPN : public testing::Test
{
    protected:
        std::vector<double> e_lda, v1_lda, v2_lda;

        void SetUp()
        {
            std::vector<double> rho  = {0.17E+01, 0.17E+01, 0.15E+01, 0.88E-01, 0.18E+04};
            std::vector<double> zeta = {0.0, 0.2, 0.5, 0.8, 1.0};

            for(int i=0;i<5;i++)
            {
                double e,v1,v2;
                XC_Functional::slater_rxc_spin(rho[i],zeta[i],e,v1,v2);
                e_lda.push_back(e);
                v1_lda.push_back(v1);
                v2_lda.push_back(v2);
            }           
        }
};

TEST_F(XCTest_SLATER_RXC_SPN, set_xc_type)
{
    std::vector<double> e_lda_ref  = {-0.880392474033,-0.888247554199,-0.892602327244,-0.378797005763,-9.99791043021};
    std::vector<double> v1_lda_ref  = {-1.17442450837,-1.24798175209,-1.28947632654,-0.532758325459,-14.0193907595};
    std::vector<double> v2_lda_ref  = {-1.17442450837,-1.09028491017,-0.894235359765,-0.2561411064,-0.688843519257};

    for (int i = 0;i<5;++i)
    {
        EXPECT_NEAR(e_lda[i],e_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v1_lda[i],v1_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v2_lda[i],v2_lda_ref[i],1.0e-8);
    }
}

class XCTest_P86_SPN : public testing::Test
{
    protected:
        std::vector<double> e_gga, v1_gga, v2_gga, v3_gga;

        void SetUp()
        {
            std::vector<double> rho  = {0.17E+01, 0.17E+01, 0.15E+01, 0.88E-01, 0.18E+04};
            std::vector<double> grho = {0.81E-11, 0.17E+01, 0.36E+02, 0.87E-01, 0.55E+00};
            std::vector<double> zeta = {0.0, 0.2, 0.5, 0.8, 1.0};

            for(int i=0;i<5;i++)
            {
                double e,v1,v2,v3;
                XC_Functional::perdew86_spin(rho[i],zeta[i],grho[i],e,v1,v2,v3);
                e_gga.push_back(e);
                v1_gga.push_back(v1);
                v2_gga.push_back(v2);
                v3_gga.push_back(v3);
            }         
        }
};

TEST_F(XCTest_P86_SPN, set_xc_type)
{
    std::vector<double> e_gga_ref  = {1.69759930415e-14,0.00308117598269,0.0408000476851,0.00290513793266,8.49954060062e-08};
    std::vector<double> v1_gga_ref  = {-1.32642497169e-14,-0.0022804667616,-0.0166964664575,-0.00685850787795,-6.30224371228e-11};
    std::vector<double> v2_gga_ref  = {-1.32642497169e-14,-0.00188528133486,-0.00317100883673,0.0160558001336,-2.36727121199e-11};
    std::vector<double> v3_gga_ref  = {0.00419160260597,0.00338159500585,0.00145591122701,0.0336591369313,3.09070721877e-07};

    for (int i = 0;i<5;++i)
    {
        EXPECT_NEAR(e_gga[i],e_gga_ref[i],1.0e-8);
        EXPECT_NEAR(v1_gga[i],v1_gga_ref[i],1.0e-8);
        EXPECT_NEAR(v2_gga[i],v2_gga_ref[i],1.0e-8);
        EXPECT_NEAR(v3_gga[i],v3_gga_ref[i],1.0e-8);
    }
}