#include "../xc_functional.h"
#include "gtest/gtest.h"
#include "../exx_info.h"

/************************************************
*  unit test of functionals
***********************************************/

// For more information of the functions, check the comment of xc_functional.h
// LDA and GGA Functionals are tested under nspin = 2

namespace ModuleBase
{
    void WARNING_QUIT(const std::string &file,const std::string &description) {exit(1);}
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

class XCTest_PBE_SPN : public testing::Test
{
    protected:
        std::vector<double> e_lda, v1_lda, v2_lda;
        std::vector<double> e_gga, v1_gga, v2_gga, v3_gga;

        void SetUp()
        {
            XC_Functional::set_xc_type("PBE");
            std::vector<double> rho  = {0.17E+01, 0.17E+01, 0.15E+01, 0.88E-01, 0.18E+04};
            std::vector<double> grho = {0.81E-11, 0.17E+01, 0.36E+02, 0.87E-01, 0.55E+00};
            std::vector<double> zeta = {0.0, 0.2, 0.5, 0.8, 1.0};

            for(int i=0;i<5;i++)
            {
                double e,ec,ex,v1,v2,v1c,v2c,v1x,v2x,v3,v4;
                XC_Functional::xc_spin(rho[i],zeta[i],e,v1,v2);
                e_lda.push_back(e);
                v1_lda.push_back(v1);
                v2_lda.push_back(v2);
                XC_Functional::gcc_spin(rho[i],zeta[i],grho[i],ec,v1c,v2c,v3);
                double r1 = rho[i] * (1+zeta[i]) / 2.0;
                double r2 = rho[i] * (1-zeta[i]) / 2.0;
                XC_Functional::gcx_spin(r1,r2,grho[i],grho[i],ex,v1x,v2x,v3,v4);
                e_gga.push_back(ec+ex);
                v1_gga.push_back(v1c+v1x);
                v2_gga.push_back(v2c+v2x);
            }
        }
};

TEST_F(XCTest_PBE_SPN, set_xc_type)
{
    EXPECT_EQ(XC_Functional::get_func_type(),2);
    std::vector<double> e_lda_ref  = {-0.9570906378,-0.9639575384,-0.961793486,-0.4179657552,-11.39185886};
    std::vector<double> v1_lda_ref  = {-1.259358955,-1.323886052,-1.352864095,-0.5692699308,-15.16995135};
    std::vector<double> v2_lda_ref  = {-1.259358955,-1.186042827,-1.011084744,-0.3739376001,-0.6598499399};
    std::vector<double> e_gga_ref  = {-5.079860763e-14,-0.0114633989,-0.1856469257,-0.005009237194,-6.500894123e-05};
    std::vector<double> v1_gga_ref = {3.971845041e-14,0.004510749212,0.07257311463,0.06143878439,-6.276046856e-15};
    std::vector<double> v2_gga_ref = {3.971845041e-14,0.01516510638,-0.1856599458,-0.1409541853,-0.09630895694};

    for (int i = 0;i<5;++i)
    {
        EXPECT_NEAR(e_lda[i],e_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v1_lda[i],v1_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v2_lda[i],v2_lda_ref[i],1.0e-8);
        EXPECT_NEAR(e_gga[i],e_gga_ref[i],1.0e-8);
        EXPECT_NEAR(v1_gga[i],v1_gga_ref[i],1.0e-8);
        EXPECT_NEAR(v2_gga[i],v2_gga_ref[i],1.0e-8);
    }
}

class XCTest_BP_SPN : public testing::Test
{
    protected:
        std::vector<double> e_lda, v1_lda, v2_lda;
        std::vector<double> e_gga, v1_gga, v2_gga, v3_gga;

        void SetUp()
        {
            XC_Functional::set_xc_type("BP");
            std::vector<double> rho  = {0.17E+01, 0.17E+01, 0.15E+01, 0.88E-01, 0.18E+04};
            std::vector<double> grho = {0.81E-11, 0.17E+01, 0.36E+02, 0.87E-01, 0.55E+00};
            std::vector<double> zeta = {0.0, 0.2, 0.5, 0.8, 1.0};

            for(int i=0;i<5;i++)
            {
                double e,ec,ex,v1,v2,v1c,v2c,v1x,v2x,v3,v4;
                XC_Functional::xc_spin(rho[i],zeta[i],e,v1,v2);
                e_lda.push_back(e);
                v1_lda.push_back(v1);
                v2_lda.push_back(v2);
                XC_Functional::gcc_spin(rho[i],zeta[i],grho[i],ec,v1c,v2c,v3);
                double r1 = rho[i] * (1+zeta[i]) / 2.0;
                double r2 = rho[i] * (1-zeta[i]) / 2.0;
                XC_Functional::gcx_spin(r1,r2,grho[i],grho[i],ex,v1x,v2x,v3,v4);
                e_gga.push_back(ec+ex);
                v1_gga.push_back(v1c+v1x);
                v2_gga.push_back(v2c+v2x);
            }
        }
};

TEST_F(XCTest_BP_SPN, set_xc_type)
{
    EXPECT_EQ(XC_Functional::get_func_type(),2);
    std::vector<double> e_lda_ref  = {-0.9565173576,-0.9631718612,-0.9599847987,-0.4169042219,-11.39281857};
    std::vector<double> v1_lda_ref  = {-1.258813137,-1.321386744,-1.348655,-0.5679440064,-15.17099378};
    std::vector<double> v2_lda_ref  = {-1.258813137,-1.18779241,-1.015533091,-0.3709293796,-0.5175834283};
    std::vector<double> e_gga_ref  = {-2.527553387e-14,-0.005791315733,-0.09840317072,-0.004480017364,-0.00633879446};
    std::vector<double> v1_gga_ref = {5.301265524e-14,0.006009233244,0.06218732775,0.04681247851,1.51250689e-11};
    std::vector<double> v2_gga_ref = {5.301265524e-14,0.01705805604,0.08693174227,-0.1388541997,-1.917431727};

    for (int i = 0;i<5;++i)
    {
        EXPECT_NEAR(e_lda[i],e_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v1_lda[i],v1_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v2_lda[i],v2_lda_ref[i],1.0e-8);
        EXPECT_NEAR(e_gga[i],e_gga_ref[i],1.0e-8);
        EXPECT_NEAR(v1_gga[i],v1_gga_ref[i],1.0e-8);
        EXPECT_NEAR(v2_gga[i],v2_gga_ref[i],1.0e-8);
    }
}

class XCTest_revPBE_SPN : public testing::Test
{
    protected:
        std::vector<double> e_lda, v1_lda, v2_lda;
        std::vector<double> e_gga, v1_gga, v2_gga, v3_gga;

        void SetUp()
        {
            XC_Functional::set_xc_type("REVPBE");
            std::vector<double> rho  = {0.17E+01, 0.17E+01, 0.15E+01, 0.88E-01, 0.18E+04};
            std::vector<double> grho = {0.81E-11, 0.17E+01, 0.36E+02, 0.87E-01, 0.55E+00};
            std::vector<double> zeta = {0.0, 0.2, 0.5, 0.8, 1.0};

            for(int i=0;i<5;i++)
            {
                double e,ec,ex,v1,v2,v1c,v2c,v1x,v2x,v3,v4;
                XC_Functional::xc_spin(rho[i],zeta[i],e,v1,v2);
                e_lda.push_back(e);
                v1_lda.push_back(v1);
                v2_lda.push_back(v2);
                XC_Functional::gcc_spin(rho[i],zeta[i],grho[i],ec,v1c,v2c,v3);
                double r1 = rho[i] * (1+zeta[i]) / 2.0;
                double r2 = rho[i] * (1-zeta[i]) / 2.0;
                XC_Functional::gcx_spin(r1,r2,grho[i],grho[i],ex,v1x,v2x,v3,v4);
                e_gga.push_back(ec+ex);
                v1_gga.push_back(v1c+v1x);
                v2_gga.push_back(v2c+v2x);
            }
        }
};

TEST_F(XCTest_revPBE_SPN, set_xc_type)
{
    EXPECT_EQ(XC_Functional::get_func_type(),2);
    std::vector<double> e_lda_ref  = {-0.9570906378,-0.9639575384,-0.961793486,-0.4179657552,-11.39185886};
    std::vector<double> v1_lda_ref  = {-1.259358955,-1.323886052,-1.352864095,-0.5692699308,-15.16995135};
    std::vector<double> v2_lda_ref  = {-1.259358955,-1.186042827,-1.011084744,-0.3739376001,-0.6598499399};
    std::vector<double> e_gga_ref  = {-5.063224319e-14,-0.01154733753,-0.2341920935,-0.006370370784,-0.0001006666667};
    std::vector<double> v1_gga_ref = {3.984893233e-14,0.004565656785,0.08385476327,0.08414965016,-7.951637783e-15};
    std::vector<double> v2_gga_ref = {3.984893233e-14,0.01556639343,-0.1143270408,-0.2468182457,-0.1491343908};

    for (int i = 0;i<5;++i)
    {
        EXPECT_NEAR(e_lda[i],e_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v1_lda[i],v1_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v2_lda[i],v2_lda_ref[i],1.0e-8);
        EXPECT_NEAR(e_gga[i],e_gga_ref[i],1.0e-8);
        EXPECT_NEAR(v1_gga[i],v1_gga_ref[i],1.0e-8);
        EXPECT_NEAR(v2_gga[i],v2_gga_ref[i],1.0e-8);
    }
}

class XCTest_PZ_SPN : public testing::Test
{
    protected:
        std::vector<double> e_lda, v1_lda, v2_lda;

        void SetUp()
        {
            XC_Functional::set_xc_type("PZ");
            std::vector<double> rho  = {0.17E+01, 0.17E+01, 0.15E+01, 0.88E-01, 0.18E+04};
            std::vector<double> zeta = {0.0, 0.2, 0.5, 0.8, 1.0};

            for(int i=0;i<5;i++)
            {
                double e,v1,v2;
                XC_Functional::xc_spin(rho[i],zeta[i],e,v1,v2);
                e_lda.push_back(e);
                v1_lda.push_back(v1);
                v2_lda.push_back(v2);
            }
        }
};

TEST_F(XCTest_PZ_SPN, set_xc_type)
{
    EXPECT_EQ(XC_Functional::get_func_type(),1);
    std::vector<double> e_lda_ref  = {-0.95651735756525513,-0.96317186117476516,-0.95998479865326902,-0.41690422190075471,-11.392818565518368};
    std::vector<double> v1_lda_ref  = {-1.2588131365004438,-1.3213867444561467,-1.3486549996482768,-0.56794400635336095,-15.170993784117668};
    std::vector<double> v2_lda_ref  = {-1.2588131365004438,-1.1877924103859072,-1.0155330910438201,-0.37092937962238148,-0.51758342826640091};

    for (int i = 0;i<5;++i)
    {
        EXPECT_NEAR(e_lda[i],e_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v1_lda[i],v1_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v2_lda[i],v2_lda_ref[i],1.0e-8);
    }
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
            std::vector<double> rho  = {-1, 0.17E+01, 0.17E+01, 0.15E+01, 0.88E-01, 0.18E+04};
            std::vector<double> zeta = {0.0, 0.0, 0.2, 0.5, 0.8, 1.0};

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
    std::vector<double> e_lda_ref  = {0.0, -0.880392474033,-0.888247554199,-0.892602327244,-0.378797005763,-9.99791043021};
    std::vector<double> v1_lda_ref  = {0.0, -1.17442450837,-1.24798175209,-1.28947632654,-0.532758325459,-14.0193907595};
    std::vector<double> v2_lda_ref  = {0.0, -1.17442450837,-1.09028491017,-0.894235359765,-0.2561411064,-0.688843519257};

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

class XCTest_PBE0_SPN : public testing::Test
{
    protected:
        std::vector<double> e_lda, v1_lda, v2_lda;
        std::vector<double> e_gga, v1_gga, v2_gga, v3_gga;

        void SetUp()
        {
            XC_Functional::set_xc_type("PBE0");
            std::vector<double> rho  = {0.17E+01, 0.17E+01, 0.15E+01, 0.88E-01, 0.18E+04};
            std::vector<double> grho = {0.81E-11, 0.17E+01, 0.36E+02, 0.87E-01, 0.55E+00};
            std::vector<double> zeta = {0.0, 0.2, 0.5, 0.8, 1.0};

            for(int i=0;i<5;i++)
            {
                double e,ec,ex,v1,v2,v1c,v2c,v1x,v2x,v3,v4;
                XC_Functional::xc_spin(rho[i],zeta[i],e,v1,v2);
                e_lda.push_back(e);
                v1_lda.push_back(v1);
                v2_lda.push_back(v2);
                XC_Functional::gcc_spin(rho[i],zeta[i],grho[i],ec,v1c,v2c,v3);
                double r1 = rho[i] * (1+zeta[i]) / 2.0;
                double r2 = rho[i] * (1-zeta[i]) / 2.0;
                XC_Functional::gcx_spin(r1,r2,grho[i],grho[i],ex,v1x,v2x,v3,v4);
                e_gga.push_back(ec+ex);
                v1_gga.push_back(v1c+v1x);
                v2_gga.push_back(v2c+v2x);
            }
        }
};

TEST_F(XCTest_PBE0_SPN, set_xc_type)
{
    EXPECT_EQ(XC_Functional::get_func_type(),4);
    std::vector<double> e_lda_ref  = {-0.7367262693,-0.7416270243,-0.7383945818,-0.323250599,-8.562036119};
    std::vector<double> v1_lda_ref  = {-0.9655397971,-1.011656606,-1.03026884,-0.4360664638,-11.39685435};
    std::vector<double> v2_lda_ref  = {-0.9655397971,-0.9132852823,-0.7874096872,-0.3099000525,-0.6598499399};
    std::vector<double> e_gga_ref  = {-3.387102535e-14,-0.007743108067,-0.1276936518,-0.003065679956,-4.873559704e-05};
    std::vector<double> v1_gga_ref = {2.647449932e-14,0.002731113185,0.0509595338,0.04805322921,-1.564140065e-11};
    std::vector<double> v2_gga_ref = {2.647449932e-14,0.01077066448,-0.13893647,-0.09075421814,-0.07223171674};

    for (int i = 0;i<5;++i)
    {
        EXPECT_NEAR(e_lda[i],e_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v1_lda[i],v1_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v2_lda[i],v2_lda_ref[i],1.0e-8);
        EXPECT_NEAR(e_gga[i],e_gga_ref[i],1.0e-8);
        EXPECT_NEAR(v1_gga[i],v1_gga_ref[i],1.0e-8);
        EXPECT_NEAR(v2_gga[i],v2_gga_ref[i],1.0e-8);
    }
}

class XCTest_PBEsol_SPN : public testing::Test
{
    protected:
        std::vector<double> e_lda, v1_lda, v2_lda;
        std::vector<double> e_gga, v1_gga, v2_gga, v3_gga;

        void SetUp()
        {
            XC_Functional::set_xc_type("PBEsol");
            std::vector<double> rho  = {0.17E+01, 0.17E+01, 0.15E+01, 0.88E-01, 0.18E+04};
            std::vector<double> grho = {0.81E-11, 0.17E+01, 0.36E+02, 0.87E-01, 0.55E+00};
            std::vector<double> zeta = {0.0, 0.2, 0.5, 0.8, 1.0};

            for(int i=0;i<5;i++)
            {
                double e,ec,ex,v1,v2,v1c,v2c,v1x,v2x,v3,v4;
                XC_Functional::xc_spin(rho[i],zeta[i],e,v1,v2);
                e_lda.push_back(e);
                v1_lda.push_back(v1);
                v2_lda.push_back(v2);
                XC_Functional::gcc_spin(rho[i],zeta[i],grho[i],ec,v1c,v2c,v3);
                double r1 = rho[i] * (1+zeta[i]) / 2.0;
                double r2 = rho[i] * (1-zeta[i]) / 2.0;
                XC_Functional::gcx_spin(r1,r2,grho[i],grho[i],ex,v1x,v2x,v3,v4);
                e_gga.push_back(ec+ex);
                v1_gga.push_back(v1c+v1x);
                v2_gga.push_back(v2c+v2x);
            }
        }
};

TEST_F(XCTest_PBEsol_SPN, set_xc_type)
{
    EXPECT_EQ(XC_Functional::get_func_type(),2);
    std::vector<double> e_lda_ref  = {-0.9570906378,-0.9639575384,-0.961793486,-0.4179657552,-11.39185886};
    std::vector<double> v1_lda_ref  = {-1.259358955,-1.323886052,-1.352864095,-0.5692699308,-15.16995135};
    std::vector<double> v2_lda_ref  = {-1.259358955,-1.186042827,-1.011084744,-0.3739376001,-0.6598499399};
    std::vector<double> e_gga_ref  = {-2.62771452e-14,-0.006047638,-0.1299218028,-0.002961686645,-6.49980662e-05};
    std::vector<double> v1_gga_ref = {2.076877742e-14,0.002195207557,0.04137297359,0.04655318731,-7.946869062e-12};
    std::vector<double> v2_gga_ref = {2.076877742e-14,0.008452420394,-0.04748419323,-0.1553803728,-0.09630827121};

    for (int i = 0;i<5;++i)
    {
        EXPECT_NEAR(e_lda[i],e_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v1_lda[i],v1_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v2_lda[i],v2_lda_ref[i],1.0e-8);
        EXPECT_NEAR(e_gga[i],e_gga_ref[i],1.0e-8);
        EXPECT_NEAR(v1_gga[i],v1_gga_ref[i],1.0e-8);
        EXPECT_NEAR(v2_gga[i],v2_gga_ref[i],1.0e-8);
    }
}

class XCTest_PBE_SPN_LibXC : public testing::Test
{
    protected:
        std::vector<double> e_gga, v1_gga, v2_gga;

        void SetUp()
        {
            XC_Functional::set_xc_type("GGA_X_PBE+GGA_C_PBE");
            std::vector<double> rho  = {0.17E+01, 0.17E+01, 0.15E+01, 0.88E-01, 0.18E+04};
            ModuleBase::Vector3<double> gdr[5];
            gdr[0] = {0.81E-11,0,0};
            gdr[1] = {0.17E+01,0,0};
            gdr[2] = {0.36E+02,0,0};
            gdr[3] = {0.87E-01,0,0};
            gdr[4] = {0.55E+00,0,0};
            std::vector<double> zeta = {0.0, 0.2, 0.5, 0.8, 1.0};

            for(int i=0;i<5;i++)
            {
                double e,v1,v2,v3,v4,v5;
                double r1 = rho[i] * (1+zeta[i]) / 2.0;
                double r2 = rho[i] * (1-zeta[i]) / 2.0;
                XC_Functional::gcxc_spin_libxc(r1,r2,gdr[i],gdr[i],e,v1,v2,v3,v4,v5);
                e_gga.push_back(e);
                v1_gga.push_back(v1+v3);
                v2_gga.push_back(v2+v4);
            }
        }
};

TEST_F(XCTest_PBE_SPN_LibXC, set_xc_type)
{
    EXPECT_EQ(XC_Functional::get_func_type(),2);
    std::vector<double> e_gga_ref  = {-1.498477706,-1.644023016,-2.248969698,-0.03699760314,-20505.34775};
    std::vector<double> v1_gga_ref = {-1.183625674,-1.326937433,-1.810935955,-0.6859883169,-15.16995245};
    std::vector<double> v2_gga_ref = {-1.183625674,-1.175091595,-1.587721737,-0.4749455756,-15.09238797};

    for (int i = 0;i<4;++i)
    {
        EXPECT_NEAR(e_gga[i],e_gga_ref[i],1.0e-8);
        EXPECT_NEAR(v1_gga[i],v1_gga_ref[i],1.0e-8);
        EXPECT_NEAR(v2_gga[i],v2_gga_ref[i],1.0e-8);
    }
}

class XCTest_PZ_SPN_LibXC : public testing::Test
{
    protected:
        std::vector<double> e_lda, v1_lda, v2_lda;

        void SetUp()
        {
            XC_Functional::set_xc_type("LDA_X+LDA_C_PZ");
            std::vector<double> rho  = {-1, 0.17E+01, 0.17E+01, 0.15E+01, 0.88E-01, 0.18E+04};
            std::vector<double> zeta = {0.0, 0.0, 0.2, 0.5, 0.8, 1.0};

            for(int i=0;i<5;i++)
            {
                double e,v1,v2;
                double r1 = rho[i] * (1+zeta[i]) / 2.0;
                double r2 = rho[i] * (1-zeta[i]) / 2.0;
                XC_Functional::xc_spin_libxc(r1,r2,e,v1,v2);
                e_lda.push_back(e);
                v1_lda.push_back(v1);
                v2_lda.push_back(v2);
            }
        }
};

TEST_F(XCTest_PZ_SPN_LibXC, set_xc_type)
{
    EXPECT_EQ(XC_Functional::get_func_type(),1);
    std::vector<double> e_lda_ref  = {0.0, -0.9565173576,-0.9631718612,-0.9599847987,-0.4169042219};
    std::vector<double> v1_lda_ref  = {0,-1.258813137,-1.321386744,-1.348655,-0.5679440064};
    std::vector<double> v2_lda_ref  = {0,-1.258813137,-1.18779241,-1.015533091,-0.3709293796};

    for (int i = 0;i<5;++i)
    {
        EXPECT_NEAR(e_lda[i],e_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v1_lda[i],v1_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v2_lda[i],v2_lda_ref[i],1.0e-8);
    }
}