#include "../xc_functional.h"
#include "gtest/gtest.h"

//a mock function of WARNING_QUIT, to avoid the uncorrected call by matrix.cpp at line 37.
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
                std::cout << std::setprecision(12) << e << " " << v << std::endl;
                XC_Functional::gcxc(rho[i],grho[i],e,v1,v2);
                e_gga.push_back(e);
                v1_gga.push_back(v1);
                v2_gga.push_back(v2);
                std::cout << std::setprecision(12) << e << " " << v1 << " " << v2 << std::endl;
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

class XCTest_OPTX : public testing::Test
{
    protected:
        std::vector<double> e_lda, v_lda;
        std::vector<double> e_gga, v1_gga, v2_gga;

        void SetUp()
        {
            std::vector<double> rho  = {0.17E+01, 0.17E+01, 0.15E+01, 0.88E-01};
            std::vector<double> grho = {0.81E-11, 0.17E+01, 0.36E+02, 0.87E-01};

            for(int i=0;i<4;i++)
            {
                double e,v1,v2;
                XC_Functional::optx(rho[i],grho[i],e,v1,v2);
                e_gga.push_back(e);
                v1_gga.push_back(v1);
                v2_gga.push_back(v2);
                std::cout << std::setprecision(12) << e << " " << v1 << " " << v2 << std::endl;
            }
        }
};

TEST_F(XCTest_OPTX, set_xc_type)
{
    std::vector<double> e_gga_ref  = { -1.5756642923000,-1.5756996817600,-1.3546578592300,-0.0358777002149};
    std::vector<double> v1_gga_ref = { -1.2358151312100,-1.2357322968800,-1.1366888917400,-0.3280597580740};
    std::vector<double> v2_gga_ref = {-4.9368029597E-15,-8.2943077227E-05,-0.002107857,-0.163514439};

    for (int i = 0;i<4;++i)
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
                std::cout << std::setprecision(12) << e << " " << v1 << " " << v2 << std::endl;
            }
        }
};

TEST_F(XCTest_WCX, set_xc_type)
{
    EXPECT_EQ(XC_Functional::get_func_type(),2);
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