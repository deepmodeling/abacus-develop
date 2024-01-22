#include "../xc_functional.h"
#include "gtest/gtest.h"
#include "../exx_info.h"

/************************************************
*  unit test of functionals
***********************************************/

// For more information of the functions, check the comment of xc_functional.h
// LDA and GGA Functionals are tested under nspin = 1

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

class XCTest_PBE : public testing::Test
{
    protected:
        std::vector<double> e_lda, v_lda;
        std::vector<double> e_gga, v1_gga, v2_gga;

        void SetUp()
        {
            XC_Functional::set_xc_type("PBE");
            std::vector<double> rho  = {0.17E+01, 0.17E+01, 0.15E+01, 0.88E-01, 0.18E+04};
            std::vector<double> grho = {0.81E-11, 0.17E+01, 0.36E+02, 0.87E-01, 0.55E+00};

            for(int i=0;i<5;i++)
            {
                double e,v,v1,v2;
                XC_Functional::xc(rho[i],e,v);
                e_lda.push_back(e);
                v_lda.push_back(v);
                XC_Functional::gcxc(rho[i],grho[i],e,v1,v2);
                e_gga.push_back(e);
                v1_gga.push_back(v1);
                v2_gga.push_back(v2);
            }
        }
};

TEST_F(XCTest_PBE, set_xc_type)
{
    EXPECT_EQ(XC_Functional::get_func_type(),2);
    std::vector<double> e_lda_ref  = {-0.9570906378,-0.9570906378,-0.9200171474,-0.3808291731,-9.124862983};
    std::vector<double> v_lda_ref  = {-1.259358955 , -1.259358955, -1.210234884,-0.4975768336,-12.12952188};
    std::vector<double> e_gga_ref  = {0.0          , -0.000103750,-0.0328708695,-0.0032277985,0.0         };
    std::vector<double> v1_gga_ref = {0.0          ,0.00021536874,0.04931694948,0.05374316118,0.0         };
    std::vector<double> v2_gga_ref = {0.0          ,-0.0002386176,-0.0025842562,-0.0825164089,0.0         };

    for (int i = 0;i<5;++i)
    {
        EXPECT_NEAR(e_lda[i],e_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v_lda[i],v_lda_ref[i],1.0e-8);
        EXPECT_NEAR(e_gga[i],e_gga_ref[i],1.0e-8);
        EXPECT_NEAR(v1_gga[i],v1_gga_ref[i],1.0e-8);
        EXPECT_NEAR(v2_gga[i],v2_gga_ref[i],1.0e-8);
    }
}

class XCTest_PBEsol : public testing::Test
{
    protected:
        std::vector<double> e_lda, v_lda;
        std::vector<double> e_gga, v1_gga, v2_gga;

        void SetUp()
        {
            XC_Functional::set_xc_type("PBEsol");
            std::vector<double> rho  = {0.17E+01, 0.17E+01, 0.15E+01, 0.88E-01, 0.18E+04};
            std::vector<double> grho = {0.81E-11, 0.17E+01, 0.36E+02, 0.87E-01, 0.55E+00};

            for(int i=0;i<5;i++)
            {
                double e,v,v1,v2;
                XC_Functional::xc(rho[i],e,v);
                e_lda.push_back(e);
                v_lda.push_back(v);
                XC_Functional::gcxc(rho[i],grho[i],e,v1,v2);
                e_gga.push_back(e);
                v1_gga.push_back(v1);
                v2_gga.push_back(v2);
            }
        }
};

TEST_F(XCTest_PBEsol, set_xc_type)
{
    EXPECT_EQ(XC_Functional::get_func_type(),2);
    std::vector<double> e_lda_ref  = {-0.9570906378,-0.9570906378,-0.9200171474,-0.3808291731,-9.124862983};
    std::vector<double> v_lda_ref  = {-1.259358955 , -1.259358955, -1.210234884,-0.4975768336,-12.12952188};
    std::vector<double> e_gga_ref  = {0.0          , 0.0003989522401,-0.008905053386,-0.001312507296,1.350759637e-08};
    std::vector<double> v1_gga_ref = {0.0          ,-0.0002447900722,0.02455798029,0.04115384844,-1.000780528e-11};
    std::vector<double> v2_gga_ref = {0.0          ,0.0004104325762,-0.001101482966,-0.04988481872,4.912388025e-08};

    for (int i = 0;i<5;++i)
    {
        EXPECT_NEAR(e_lda[i],e_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v_lda[i],v_lda_ref[i],1.0e-8);
        EXPECT_NEAR(e_gga[i],e_gga_ref[i],1.0e-8);
        EXPECT_NEAR(v1_gga[i],v1_gga_ref[i],1.0e-8);
        EXPECT_NEAR(v2_gga[i],v2_gga_ref[i],1.0e-8);
    }
}

class XCTest_BP : public testing::Test
{
    protected:
        std::vector<double> e_lda, v_lda;
        std::vector<double> e_gga, v1_gga, v2_gga;

        void SetUp()
        {
            XC_Functional::set_xc_type("BP");
            std::vector<double> rho  = {0.17E+01, 0.17E+01, 0.15E+01, 0.88E-01, 0.18E+04};
            std::vector<double> grho = {0.81E-11, 0.17E+01, 0.36E+02, 0.87E-01, 0.55E+00};

            for(int i=0;i<5;i++)
            {
                double e,v,v1,v2;
                XC_Functional::xc(rho[i],e,v);
                e_lda.push_back(e);
                v_lda.push_back(v);
                XC_Functional::gcxc(rho[i],grho[i],e,v1,v2);
                e_gga.push_back(e);
                v1_gga.push_back(v1);
                v2_gga.push_back(v2);
            }
        }
};

TEST_F(XCTest_BP, set_xc_type)
{
    EXPECT_EQ(XC_Functional::get_func_type(),2);
    std::vector<double> e_lda_ref  = {-0.9565173576,-0.9565173576,-0.9194417358,-0.38104767571,-9.125575429};
    std::vector<double> v_lda_ref  = {-1.2588131365,-1.2588131365,-1.2096658201,-0.49757720968,-12.13038680};
    std::vector<double> e_gga_ref  = {0.0          ,-0.0012525632,-0.0457062154,-0.0035030371,0.0         };
    std::vector<double> v1_gga_ref = {0.0          ,0.00118557755,0.04272756138,0.04137909619,0.0         };
    std::vector<double> v2_gga_ref = {0.0,-0.0016488916,-0.0027047176,-0.0753933855,0.0         };
    for (int i = 0;i<5;++i)
    {
        EXPECT_NEAR(e_lda[i],e_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v_lda[i],v_lda_ref[i],1.0e-8);
        EXPECT_NEAR(e_gga[i],e_gga_ref[i],1.0e-7);
        EXPECT_NEAR(v1_gga[i],v1_gga_ref[i],1.0e-7);
        EXPECT_NEAR(v2_gga[i],v2_gga_ref[i],1.0e-7);
    }
}

class XCTest_revPBE : public testing::Test
{
    protected:
        std::vector<double> e_lda, v_lda;
        std::vector<double> e_gga, v1_gga, v2_gga;

        void SetUp()
        {
            XC_Functional::set_xc_type("REVPBE");
            std::vector<double> rho  = {0.17E+01, 0.17E+01, 0.15E+01, 0.88E-01, 0.18E+04};
            std::vector<double> grho = {0.81E-11, 0.17E+01, 0.36E+02, 0.87E-01, 0.55E+00};

            for(int i=0;i<5;i++)
            {
                double e,v,v1,v2;
                XC_Functional::xc(rho[i],e,v);
                e_lda.push_back(e);
                v_lda.push_back(v);
                XC_Functional::gcxc(rho[i],grho[i],e,v1,v2);
                e_gga.push_back(e);
                v1_gga.push_back(v1);
                v2_gga.push_back(v2);
            }
        }
};

TEST_F(XCTest_revPBE, set_xc_type)
{
    EXPECT_EQ(XC_Functional::get_func_type(),2);
    std::vector<double> e_lda_ref  = {-0.957090637829,-0.957090637829,-0.920017147403,-0.38082917306,-9.12486298288};
    std::vector<double> v_lda_ref  = {-1.25935895511,-1.25935895511,-1.21023488409,-0.497576833567,-12.1295218848};
    std::vector<double> e_gga_ref  = {0,-0.000107434044223,-0.0352560914762,-0.00398961883151,-1.41140669279e-12};
    std::vector<double> v1_gga_ref = {0,0.000224008875928,0.055112094687,0.0769311998392,-1.04475574595e-15};
    std::vector<double> v2_gga_ref = {0,-0.000247264687907,-0.00283161066613,-0.108863890698,-7.3118287463e-16};

    for (int i = 0;i<5;++i)
    {
        EXPECT_NEAR(e_lda[i],e_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v_lda[i],v_lda_ref[i],1.0e-8);
        EXPECT_NEAR(e_gga[i],e_gga_ref[i],1.0e-8);
        EXPECT_NEAR(v1_gga[i],v1_gga_ref[i],1.0e-8);
        EXPECT_NEAR(v2_gga[i],v2_gga_ref[i],1.0e-8);
    }
}

class XCTest_WC : public testing::Test
{
    protected:
        std::vector<double> e_lda, v_lda;
        std::vector<double> e_gga, v1_gga, v2_gga;

        void SetUp()
        {
            XC_Functional::set_xc_type("WC");
            std::vector<double> rho  = {0.17E+01, 0.17E+01, 0.15E+01, 0.88E-01, 0.18E+04};
            std::vector<double> grho = {0.81E-11, 0.17E+01, 0.36E+02, 0.87E-01, 0.55E+00};

            for(int i=0;i<5;i++)
            {
                double e,v,v1,v2;
                XC_Functional::xc(rho[i],e,v);
                e_lda.push_back(e);
                v_lda.push_back(v);
                XC_Functional::gcxc(rho[i],grho[i],e,v1,v2);
                e_gga.push_back(e);
                v1_gga.push_back(v1);
                v2_gga.push_back(v2);
            }
        }
};

TEST_F(XCTest_WC, set_xc_type)
{
    EXPECT_EQ(XC_Functional::get_func_type(),2);
    std::vector<double> e_lda_ref  = {-0.957090637829,-0.957090637829,-0.920017147403,-0.38082917306,-9.12486298288};
    std::vector<double> v_lda_ref  = {-1.25935895511,-1.25935895511,-1.21023488409,-0.497576833567,-12.1295218848};
    std::vector<double> e_gga_ref  = {0,-8.85638848139e-05,-0.0246877948944,-0.00174638618626,-1.41140669277e-12};
    std::vector<double> v1_gga_ref = {0,0.000179915620782,0.03203380922,0.0380174720271,-1.04476045089e-15};
    std::vector<double> v2_gga_ref = {0,-0.000203094652778,-0.00181685045134,-0.053558848389,-7.19634374487e-16};

    for (int i = 0;i<5;++i)
    {
        EXPECT_NEAR(e_lda[i],e_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v_lda[i],v_lda_ref[i],1.0e-8);
        EXPECT_NEAR(e_gga[i],e_gga_ref[i],1.0e-8);
        EXPECT_NEAR(v1_gga[i],v1_gga_ref[i],1.0e-8);
        EXPECT_NEAR(v2_gga[i],v2_gga_ref[i],1.0e-8);
    }
}

class XCTest_BLYP : public testing::Test
{
    protected:
        std::vector<double> e_lda, v_lda;
        std::vector<double> e_gga, v1_gga, v2_gga;

        void SetUp()
        {
            XC_Functional::set_xc_type("BLYP");
            std::vector<double> rho  = {0.17E+01, 0.17E+01, 0.15E+01, 0.88E-01, 0.18E+04};
            std::vector<double> grho = {0.81E-11, 0.17E+01, 0.36E+02, 0.87E-01, 0.55E+00};

            for(int i=0;i<5;i++)
            {
                double e,v,v1,v2;
                XC_Functional::xc(rho[i],e,v);
                e_lda.push_back(e);
                v_lda.push_back(v);
                XC_Functional::gcxc(rho[i],grho[i],e,v1,v2);
                e_gga.push_back(e);
                v1_gga.push_back(v1);
                v2_gga.push_back(v2);
            }
        }
};

TEST_F(XCTest_BLYP, set_xc_type)
{
    EXPECT_EQ(XC_Functional::get_func_type(),2);
    std::vector<double> e_lda_ref  = {-0.931174253662,-0.931174253662,-0.894576862875,-0.361975213747,-9.04968167694};
    std::vector<double> v_lda_ref  = {-1.22956821264,-1.22956821264,-1.18106120453,-0.477504409374,-12.0451233615};
    std::vector<double> e_gga_ref  = {0,-0.00412741708887,-0.0830008150996,-0.0054607659128,-1.3232827772e-07};
    std::vector<double> v1_gga_ref = {0,0.00309536590829,0.0499783203384,0.021894292547,9.79056031428e-11};
    std::vector<double> v2_gga_ref = {0,-0.00478506213734,-0.00391088520577,-0.0816833688647,-4.8119373714e-07};

    for (int i = 0;i<5;++i)
    {
        EXPECT_NEAR(e_lda[i],e_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v_lda[i],v_lda_ref[i],1.0e-8);
        EXPECT_NEAR(e_gga[i],e_gga_ref[i],1.0e-8);
        EXPECT_NEAR(v1_gga[i],v1_gga_ref[i],1.0e-8);
        EXPECT_NEAR(v2_gga[i],v2_gga_ref[i],1.0e-8);
    }
}

class XCTest_PW91 : public testing::Test
{
    protected:
        std::vector<double> e_lda, v_lda;
        std::vector<double> e_gga, v1_gga, v2_gga;

        void SetUp()
        {
            XC_Functional::set_xc_type("PW91");
            std::vector<double> rho  = {0.17E+01, 0.17E+01, 0.15E+01, 0.88E-01, 0.18E+04};
            std::vector<double> grho = {0.81E-11, 0.17E+01, 0.36E+02, 0.87E-01, 0.55E+00};

            for(int i=0;i<5;i++)
            {
                double e,v,v1,v2;
                XC_Functional::xc(rho[i],e,v);
                e_lda.push_back(e);
                v_lda.push_back(v);
                XC_Functional::gcxc(rho[i],grho[i],e,v1,v2);
                e_gga.push_back(e);
                v1_gga.push_back(v1);
                v2_gga.push_back(v2);
            }
        }
};

TEST_F(XCTest_PW91, set_xc_type)
{
    EXPECT_EQ(XC_Functional::get_func_type(),2);
    std::vector<double> e_lda_ref  = {-0.957090637829,-0.957090637829,-0.920017147403,-0.38082917306,-9.12486298288};
    std::vector<double> v_lda_ref  = {-1.25935895511 , -1.25935895511, -1.21023488409,-0.497576833567,-12.1295218848};
    std::vector<double> e_gga_ref  = {0.0          , 6.43004695777e-05,-0.0438528949486,-0.00342032994996,6.51882950017e-08};
    std::vector<double> v1_gga_ref = {0.0          ,0.00179434688187,0.0464532357723,0.0467404236193,0.0         };
    std::vector<double> v2_gga_ref = {0.0          ,-0.00132500065922,-0.00280082132252,-0.0793620638752,2.37048336664e-07};

    for (int i = 0;i<5;++i)
    {
        EXPECT_NEAR(e_lda[i],e_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v_lda[i],v_lda_ref[i],1.0e-8);
        EXPECT_NEAR(e_gga[i],e_gga_ref[i],1.0e-8);
        EXPECT_NEAR(v1_gga[i],v1_gga_ref[i],1.0e-8);
        EXPECT_NEAR(v2_gga[i],v2_gga_ref[i],1.0e-8);
    }
}

class XCTest_OLYP : public testing::Test
{
    protected:
        std::vector<double> e_lda, v_lda;
        std::vector<double> e_gga, v1_gga, v2_gga;

        void SetUp()
        {
            XC_Functional::set_xc_type("OLYP");
            std::vector<double> rho  = {0.17E+01, 0.17E+01, 0.15E+01, 0.88E-01};
            std::vector<double> grho = {0.81E-11, 0.17E+01, 0.36E+02, 0.87E-01};

            for(int i=0;i<4;i++)
            {
                double e,v,v1,v2;
                XC_Functional::xc(rho[i],e,v);
                e_lda.push_back(e);
                v_lda.push_back(v);
                XC_Functional::gcxc(rho[i],grho[i],e,v1,v2);
                e_gga.push_back(e);
                v1_gga.push_back(v1);
                v2_gga.push_back(v2);
            }
        }
};

TEST_F(XCTest_OLYP, set_xc_type)
{
    EXPECT_EQ(XC_Functional::get_func_type(),2);
    std::vector<double> e_lda_ref  = {-0.0497167797343,-0.0497167797343,-0.0491381240345,-0.0334672153176};
    std::vector<double> v_lda_ref  = {-0.0542915807395,-0.0542915807395,-0.0538095527425,-0.0394937448009};
    std::vector<double> e_gga_ref  = {0,-1.57545926183,-1.34838480556,-0.0344400770416};
    std::vector<double> v1_gga_ref = {0,-1.2359683886,-1.14365644338,-0.352882289038};
    std::vector<double> v2_gga_ref = {0,0.000199903898453,-0.00175935413007,-0.130465630292};

    for (int i = 0;i<4;++i)
    {
        EXPECT_NEAR(e_lda[i],e_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v_lda[i],v_lda_ref[i],1.0e-8);
        EXPECT_NEAR(e_gga[i],e_gga_ref[i],1.0e-8);
        EXPECT_NEAR(v1_gga[i],v1_gga_ref[i],1.0e-8);
        EXPECT_NEAR(v2_gga[i],v2_gga_ref[i],1.0e-8);
    }
}


class XCTest_HCTH : public testing::Test
{
    protected:
        std::vector<double> e_lda, v_lda;
        std::vector<double> e_gga, v1_gga, v2_gga;

        void SetUp()
        {
            XC_Functional::set_xc_type("HCTH");
            std::vector<double> rho  = {0.17E+01, 0.17E+01, 0.15E+01, 0.88E-01};
            std::vector<double> grho = {0.81E-11, 0.17E+01, 0.36E+02, 0.87E-01};

            for(int i=0;i<4;i++)
            {
                double e,v,v1,v2;
                XC_Functional::xc(rho[i],e,v);
                e_lda.push_back(e);
                v_lda.push_back(v);
                XC_Functional::gcxc(rho[i],grho[i],e,v1,v2);
                e_gga.push_back(e);
                v1_gga.push_back(v1);
                v2_gga.push_back(v2);
            }
        }
};

TEST_F(XCTest_HCTH, set_xc_type)
{
    EXPECT_EQ(XC_Functional::get_func_type(),2);
    std::vector<double> e_lda_ref  = {0,0,0,0};
    std::vector<double> v_lda_ref  = {0,0,0,0};
    std::vector<double> e_gga_ref  = {0,-1.69759307165,-1.43450950992,-0.0379281863915};
    std::vector<double> v1_gga_ref = {0,-1.32696036531,-1.22452673027,-0.417946684263};
    std::vector<double> v2_gga_ref = {0,0.00296525371089,-0.00118057767393,-0.112470272644};

    for (int i = 0;i<4;++i)
    {
        EXPECT_NEAR(e_lda[i],e_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v_lda[i],v_lda_ref[i],1.0e-8);
        EXPECT_NEAR(e_gga[i],e_gga_ref[i],1.0e-8);
        EXPECT_NEAR(v1_gga[i],v1_gga_ref[i],1.0e-8);
        EXPECT_NEAR(v2_gga[i],v2_gga_ref[i],1.0e-8);
    }
}

class XCTest_PWLDA : public testing::Test
{
    protected:
        std::vector<double> e_lda, v_lda;

        void SetUp()
        {
            XC_Functional::set_xc_type("PWLDA");
            std::vector<double> rho  = {0.17E+01, 0.17E+01, 0.15E+01, 0.88E-01, 0.18E+04};

            for(int i=0;i<5;i++)
            {
                double e,v;
                XC_Functional::xc(rho[i],e,v);
                e_lda.push_back(e);
                v_lda.push_back(v);
            }
        }
};

TEST_F(XCTest_PWLDA, set_xc_type)
{
    EXPECT_EQ(XC_Functional::get_func_type(),1);
    std::vector<double> e_lda_ref  = {-0.9570906378,-0.9570906378,-0.9200171474,-0.3808291731,-9.124862983};
    std::vector<double> v_lda_ref  = {-1.259358955 , -1.259358955, -1.210234884,-0.4975768336,-12.12952188};

    for (int i = 0;i<5;++i)
    {
        EXPECT_NEAR(e_lda[i],e_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v_lda[i],v_lda_ref[i],1.0e-8);
    }
}

class XCTest_PZ : public testing::Test
{
    protected:
        std::vector<double> e_lda, v_lda;

        void SetUp()
        {
            XC_Functional::set_xc_type("PZ");
            std::vector<double> rho  = {0.17E+01, 0.17E+01, 0.15E+01, 0.88E-01, 0.18E+04};

            for(int i=0;i<5;i++)
            {
                double e,v;
                XC_Functional::xc(rho[i],e,v);
                e_lda.push_back(e);
                v_lda.push_back(v);
            }
        }
};

TEST_F(XCTest_PZ, set_xc_type)
{
    EXPECT_EQ(XC_Functional::get_func_type(),1);
    std::vector<double> e_lda_ref  = {-0.956517357565,-0.956517357565,-0.91944173584,-0.381047675708,-9.1255754288};
    std::vector<double> v_lda_ref  = {-1.258813137,-1.258813137,-1.20966582,-0.49757721,-12.1303868};

    for (int i = 0;i<5;++i)
    {
        EXPECT_NEAR(e_lda[i],e_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v_lda[i],v_lda_ref[i],1.0e-8);
    }
}

class XCTest_SLATER1 : public testing::Test
{
    protected:
        std::vector<double> e_lda, v_lda;

        void SetUp()
        {
            std::vector<double> rho  = {0.17E+01, 0.17E+01, 0.15E+01, 0.88E-01, 0.18E+04};

            for(int i=0;i<5;i++)
            {
                double e,v;
                XC_Functional::slater1(rho[i],e,v);
                e_lda.push_back(e);
                v_lda.push_back(v);
            }
        }
};

TEST_F(XCTest_SLATER1, set_xc_type)
{
    std::vector<double> e_lda_ref  = {-0.404263494,-0.404263494,-0.458165293,-7.809635681,-0.000381804};
    std::vector<double> v_lda_ref  = {-0.539017992,-0.539017992,-0.610887058,-10.41284757,-0.000509073};

    for (int i = 0;i<5;++i)
    {
        EXPECT_NEAR(e_lda[i],e_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v_lda[i],v_lda_ref[i],1.0e-8);
    }
}

class XCTest_SLATER_RXC : public testing::Test
{
    protected:
        std::vector<double> e_lda, v_lda;

        void SetUp()
        {
            std::vector<double> rho  = {0.17E+01, 0.17E+01, 0.15E+01, 0.88E-01, 0.18E+04};

            for(int i=0;i<5;i++)
            {
                double e,v;
                XC_Functional::slater_rxc(rho[i],e,v);
                e_lda.push_back(e);
                v_lda.push_back(v);
            }
        }
};

TEST_F(XCTest_SLATER_RXC, set_xc_type)
{
    std::vector<double> e_lda_ref  = {-0.269496811,-0.269496811,-0.305425791,-5.1198857,-0.000254536};
    std::vector<double> v_lda_ref  = {-0.359320959,-0.359320959,-0.407222564,-6.769681567,-0.000339382};

    for (int i = 0;i<5;++i)
    {
        EXPECT_NEAR(e_lda[i],e_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v_lda[i],v_lda_ref[i],1.0e-8);
    }
}

class XCTest_PW : public testing::Test
{
    protected:
        std::vector<double> e_lda, v_lda;

        void SetUp()
        {
            std::vector<double> rho  = {0.1, 0.2, 101, 110};

            for(int i=0;i<4;i++)
            {
                double e,v;
                XC_Functional::pw(rho[i],1,e,v);
                e_lda.push_back(e);
                v_lda.push_back(v);
            }
        }
};

TEST_F(XCTest_PW, set_xc_type)
{
    std::vector<double> e_lda_ref  = {-0.1208055896280259,-0.1009063676832011,-0.0028726241745994537,-0.0026920464558317734};
    std::vector<double> v_lda_ref  = {-0.13053328412744325,-0.11030492316729694,-0.0035935897272457141,-0.0033812515019294776};

    for (int i = 0;i<4;++i)
    {
        EXPECT_NEAR(e_lda[i],e_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v_lda[i],v_lda_ref[i],1.0e-8);
    }
}


class XCTest_LYP : public testing::Test
{
    protected:
        std::vector<double> e_lda, v_lda;

        void SetUp()
        {
            std::vector<double> rho  = {0.17E+01, 0.17E+01, 0.15E+01, 0.88E-01, 0.18E+04};

            for(int i=0;i<5;i++)
            {
                double e,v;
                XC_Functional::lyp(rho[i],e,v);
                e_lda.push_back(e);
                v_lda.push_back(v);
            }
        }
};

TEST_F(XCTest_LYP, set_xc_type)
{
    std::vector<double> e_lda_ref  = {-0.02989705,-0.02989705,-0.032151081,-0.063993353,-4.8517572814898958e-05};
    std::vector<double> v_lda_ref  = {-0.035869963,-0.035869963,-0.038174479,-0.065204797,-6.4674142397147037e-05};

    for (int i = 0;i<5;++i)
    {
        EXPECT_NEAR(e_lda[i],e_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v_lda[i],v_lda_ref[i],1.0e-8);
    }
}

class XCTest_VWN : public testing::Test
{
    protected:
        std::vector<double> e_lda, v_lda;

        void SetUp()
        {
            std::vector<double> rho  = {0.17E+01, 0.17E+01, 0.15E+01, 0.88E-01, 0.18E+04};

            for(int i=0;i<5;i++)
            {
                double e,v;
                XC_Functional::vwn(rho[i],e,v);
                e_lda.push_back(e);
                v_lda.push_back(v);
            }
        }
};

TEST_F(XCTest_VWN, set_xc_type)
{
    std::vector<double> e_lda_ref  = { -0.0481692845703,-0.0481692845703,-0.0508581491890,-0.1250586358990,-0.0002168020098};
    std::vector<double> v_lda_ref  = { -0.0552381037021,-0.0552381037021,-0.0581101146837,-0.1347585962230,-0.0002868652922};

    for (int i = 0;i<5;++i)
    {
        EXPECT_NEAR(e_lda[i],e_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v_lda[i],v_lda_ref[i],1.0e-8);
    }
}

class XCTest_WIGNER : public testing::Test
{
    protected:
        std::vector<double> e_lda, v_lda;

        void SetUp()
        {
            std::vector<double> rho  = {0.17E+01, 0.17E+01, 0.15E+01, 0.88E-01, 0.18E+04};

            for(int i=0;i<5;i++)
            {
                double e,v;
                XC_Functional::wigner(rho[i],e,v);
                e_lda.push_back(e);
                v_lda.push_back(v);
            }
        }
};

TEST_F(XCTest_WIGNER, set_xc_type)
{
    std::vector<double> e_lda_ref  = {-0.046226266609,-0.046226266609,-0.047220614485,-0.055675743921,-0.000242863497};
    std::vector<double> v_lda_ref  = {-0.048984187643,-0.048984187643,-0.049759880946,-0.055882739698,-0.000323468803};

    for (int i = 0;i<5;++i)
    {
        EXPECT_NEAR(e_lda[i],e_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v_lda[i],v_lda_ref[i],1.0e-8);
    }
}

class XCTest_HL : public testing::Test
{
    protected:
        std::vector<double> e_lda, v_lda;

        void SetUp()
        {
            std::vector<double> rho  = {0.17E+01, 0.17E+01, 0.15E+01, 0.88E-01, 0.18E+04};

            for(int i=0;i<5;i++)
            {
                double e,v;
                XC_Functional::hl(rho[i],e,v);
                e_lda.push_back(e);
                v_lda.push_back(v);
            }
        }
};

TEST_F(XCTest_HL, set_xc_type)
{
    std::vector<double> e_lda_ref  = {-0.0516082763238,-0.0516082763238,-0.0541421102571,-0.1158270158400,-0.0001959621612};
    std::vector<double> v_lda_ref  = {-0.0583140751521,-0.0583140751521,-0.0609311295248,-0.1232802590140,-0.0002609805565 };

    for (int i = 0;i<5;++i)
    {
        EXPECT_NEAR(e_lda[i],e_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v_lda[i],v_lda_ref[i],1.0e-8);
    }
}

class XCTest_GL : public testing::Test
{
    protected:
        std::vector<double> e_lda, v_lda;

        void SetUp()
        {
            std::vector<double> rho  = {0.17E+01, 0.17E+01, 0.15E+01, 0.88E-01, 0.18E+04};

            for(int i=0;i<5;i++)
            {
                double e,v;
                XC_Functional::gl(rho[i],e,v);
                e_lda.push_back(e);
                v_lda.push_back(v);
            }
        }
};

TEST_F(XCTest_GL, set_xc_type)
{
    std::vector<double> e_lda_ref  = {-0.0588659395493,-0.0588659395493,-0.0623311764511,-0.1512549418480,-0.0001577756901};
    std::vector<double> v_lda_ref  = {-0.0679980665055,-0.0679980665055,-0.0716536813685,-0.1622283251780,-0.0002102349565};

    for (int i = 0;i<5;++i)
    {
        EXPECT_NEAR(e_lda[i],e_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v_lda[i],v_lda_ref[i],1.0e-8);
    }
}

class XCTest_OPTX : public testing::Test
{
    protected:
        std::vector<double> e_lda, v_lda;
        std::vector<double> e_gga, v1_gga, v2_gga;

        void SetUp()
        {
            std::vector<double> rho  = {0.0, 0.17E+01, 0.17E+01, 0.15E+01, 0.88E-01};
            std::vector<double> grho = {0.0, 0.81E-11, 0.17E+01, 0.36E+02, 0.87E-01};

            for(int i=0;i<5;i++)
            {
                double e,v1,v2;
                XC_Functional::optx(rho[i],grho[i],e,v1,v2);
                e_gga.push_back(e);
                v1_gga.push_back(v1);
                v2_gga.push_back(v2);
            }
        }
};

TEST_F(XCTest_OPTX, set_xc_type)
{
    std::vector<double> e_gga_ref  = { 0.0, -1.5756642923000,-1.5756996817600,-1.3546578592300,-0.0358777002149};
    std::vector<double> v1_gga_ref = { 0.0, -1.2358151312100,-1.2357322968800,-1.1366888917400,-0.3280597580740};
    std::vector<double> v2_gga_ref = {-0.0, 4.9368029597E-15,-8.2943077227E-05,-0.002107857,-0.163514439};

    for (int i = 0;i<5;++i)
    {
        EXPECT_NEAR(e_gga[i],e_gga_ref[i],1.0e-8);
        EXPECT_NEAR(v1_gga[i],v1_gga_ref[i],1.0e-8);
        EXPECT_NEAR(v2_gga[i],v2_gga_ref[i],1.0e-8);
    }
}

class XCTest_WCX : public testing::Test
{
    protected:
        std::vector<double> e_lda, v_lda;
        std::vector<double> e_gga, v1_gga, v2_gga;

        void SetUp()
        {
            std::vector<double> rho  = {0.17E+01, 0.17E+01, 0.15E+01, 0.88E-01, 0.18E+04};
            std::vector<double> grho = {0.81E-11, 0.17E+01, 0.36E+02, 0.87E-01, 0.55E+00};

            for(int i=0;i<5;i++)
            {
                double e,v1,v2;
                XC_Functional::wcx(rho[i],grho[i],e,v1,v2);
                e_gga.push_back(e);
                v1_gga.push_back(v1);
                v2_gga.push_back(v2);
            }
        }
};

TEST_F(XCTest_WCX, set_xc_type)
{
    std::vector<double> e_gga_ref  = {0.0          ,-0.0035227366077,-0.0734928155169,-0.0052181684189,-0.0000001063768 };
    std::vector<double> v1_gga_ref = {0.0          ,0.0027230912249,0.0436861496717,0.0272641442833,0.0000000000788 };
    std::vector<double> v2_gga_ref = {-4.1745209791E-03,-4.1145164232E-03,-3.4066592749E-03,-8.0662091283E-02,-3.8681965771E-07};

    for (int i = 0;i<5;++i)
    {
        EXPECT_NEAR(e_gga[i],e_gga_ref[i],1.0e-8);
        EXPECT_NEAR(v1_gga[i],v1_gga_ref[i],1.0e-8);
        EXPECT_NEAR(v2_gga[i],v2_gga_ref[i],1.0e-8);
    }
}

class XCTest_PBE0 : public testing::Test
{
    protected:
        std::vector<double> e_lda, v_lda;
        std::vector<double> e_gga, v1_gga, v2_gga;

        void SetUp()
        {
            XC_Functional::set_xc_type("PBE0");
            XC_Functional::get_hybrid_alpha(0.5);
            std::vector<double> rho  = {0.17E+01, 0.17E+01, 0.15E+01, 0.88E-01, 0.18E+04};
            std::vector<double> grho = {0.81E-11, 0.17E+01, 0.36E+02, 0.87E-01, 0.55E+00};

            for(int i=0;i<5;i++)
            {
                double e,v,v1,v2;
                XC_Functional::xc(rho[i],e,v);
                e_lda.push_back(e);
                v_lda.push_back(v);
                XC_Functional::gcxc(rho[i],grho[i],e,v1,v2);
                e_gga.push_back(e);
                v1_gga.push_back(v1);
                v2_gga.push_back(v2);
            }
        }
};

TEST_F(XCTest_PBE0, set_xc_type)
{
    EXPECT_EQ(XC_Functional::get_func_type(),4);
    std::vector<double> e_lda_ref  = {-0.516361900865,-0.516361900865,-0.497297777983,-0.216575173845,-4.63279938014};
    std::vector<double> v_lda_ref  = {-0.671720639159,-0.671720639159,-0.646609058192,-0.27857150128,-6.14010374778};
    std::vector<double> e_gga_ref  = {0,0.00166521130532,0.00796707553969,0.00012199188186,5.31869959642e-08};
    std::vector<double> v1_gga_ref = {0,-0.00116390343176,0.0188323045143,0.0322482444605,-3.93988214247e-11};
    std::vector<double> v2_gga_ref = {0,0.00183640210141,-0.000497223695012,-0.0277065830152,1.9340982813e-07};

    for (int i = 0;i<5;++i)
    {
        EXPECT_NEAR(e_lda[i],e_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v_lda[i],v_lda_ref[i],1.0e-8);
        EXPECT_NEAR(e_gga[i],e_gga_ref[i],1.0e-8);
        EXPECT_NEAR(v1_gga[i],v1_gga_ref[i],1.0e-8);
        EXPECT_NEAR(v2_gga[i],v2_gga_ref[i],1.0e-8);
    }
}

class XCTest_PBE_LibXC : public testing::Test
{
    protected:
        std::vector<double> e_lda, v_lda;
        std::vector<double> e_gga, v1_gga, v2_gga;

        void SetUp()
        {
            XC_Functional::set_xc_type("GGA_X_PBE+GGA_C_PBE");
            std::vector<double> rho  = {0.17E+01, 0.17E+01, 0.15E+01, 0.88E-01, 0.18E+04};
            std::vector<double> grho = {0.81E-11, 0.17E+01, 0.36E+02, 0.87E-01, 0.55E+00};

            for(int i=0;i<5;i++)
            {
                double e,v,v1,v2;
                XC_Functional::xc(rho[i],e,v);
                e_lda.push_back(e);
                v_lda.push_back(v);
                XC_Functional::gcxc_libxc(rho[i],grho[i],e,v1,v2);
                e_gga.push_back(e);
                v1_gga.push_back(v1);
                v2_gga.push_back(v2);
            }
        }
};

TEST_F(XCTest_PBE_LibXC, set_xc_type)
{
    EXPECT_EQ(XC_Functional::get_func_type(),2);
    std::vector<double> e_lda_ref  = {-0.9570906378,-0.9570906378,-0.9200171474,-0.3808291731,-9.124862983};
    std::vector<double> v_lda_ref  = {-1.259358955 , -1.259358955, -1.210234884,-0.4975768336,-12.12952188};
    std::vector<double> e_gga_ref  = {0.0          , -1.6271573478357044,-1.412896188863435,-0.036740758114697646,0.0         };
    std::vector<double> v1_gga_ref = {0.0          ,-1.2591432392430806,-1.1609175516612873,-0.44383341607833338,0.0         };
    std::vector<double> v2_gga_ref = {0.0          ,-0.0002386176,-0.0025842562,-0.082516520555534323,0.0         };

    for (int i = 0;i<4;++i)
    {
        EXPECT_NEAR(e_lda[i],e_lda_ref[i],1.0e-8);
        EXPECT_NEAR(v_lda[i],v_lda_ref[i],1.0e-8);
        EXPECT_NEAR(e_gga[i],e_gga_ref[i],1.0e-8);
        EXPECT_NEAR(v1_gga[i],v1_gga_ref[i],1.0e-8);
        EXPECT_NEAR(v2_gga[i],v2_gga_ref[i],1.0e-8);
    }
}

/*
//for printing results
            std::cout << std::setprecision(10);
            for(int i=0;i<5;i++)
            {
                std::cout << e_lda[i] << ",";
            }
            std::cout << std::endl;
*/