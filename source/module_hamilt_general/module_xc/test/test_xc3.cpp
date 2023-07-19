#include "../xc_functional.h"
#include "gtest/gtest.h"
#include "../exx_info.h"
#include "xc3_mock.h"
#include "module_base/matrix.h"

/************************************************
*  unit test of functionals
***********************************************/

// For more information of the functions, check the comment of xc_functional.h
// Two functions are tested:
// gradcorr, which calculates the gradient part of GGA functional
// gradwfc, which is used to obtain the derivative of wavefunction

class XCTest_GRADCORR : public testing::Test
{
    protected:

        double et1 = 0, vt1 = 0;
        ModuleBase::matrix v1;
        std::vector<double> stress1;

        double et2 = 0, vt2 = 0;
        ModuleBase::matrix v2;
        std::vector<double> stress2;

        double et4 = 0, vt4 = 0;
        ModuleBase::matrix v4;
        std::vector<double> stress4;

        void SetUp()
        {
            ModulePW::PW_Basis rhopw;
            UnitCell ucell;
            Charge chr;

            rhopw.nrxx = 5;
            rhopw.npw = 5;
            rhopw.nmaxgr = 5;
            rhopw.gcar = new ModuleBase::Vector3<double> [5];

            ucell.tpiba = 1;
            ucell.magnet.lsign_ = true;
            ucell.cal_ux();

            chr.rho = new double*[4];
            chr.rho[0] = new double[5];
            chr.rho[1] = new double[5];
            chr.rho[2] = new double[5];
            chr.rho[3] = new double[5];
            chr.rhog = new std::complex<double>*[2];
            chr.rhog[0] = new std::complex<double>[5];
            chr.rhog[1] = new std::complex<double>[5];

            chr.rho_core = new double[5];
            chr.rhog_core = new std::complex<double>[5];

            for(int i=0;i<5;i++)
            {
                chr.rho[0][i] = double(i);
                chr.rho[1][i] = 0.1*double(i);
                chr.rho[2][i] = chr.rho[0][i];
                chr.rho[3][i] = chr.rho[1][i];
                chr.rhog[0][i] = chr.rho[0][i];
                chr.rhog[1][i] = chr.rho[1][i];
                chr.rho_core[i] = 0;
                chr.rhog_core[i] = 0;
                rhopw.gcar[i]= 1;
            }

            v1.create(1,5);
            v1.zero_out();
            v2.create(2,5);
            v2.zero_out();
            v4.create(4,5);
            v4.zero_out();

            XC_Functional::set_xc_type("PBE");

            GlobalV::NSPIN = 1;
            XC_Functional::gradcorr(et1,vt1,v1,&chr,&rhopw,&ucell,stress1,false);
            XC_Functional::gradcorr(et1,vt1,v1,&chr,&rhopw,&ucell,stress1,true);

            GlobalV::NSPIN = 2;
            XC_Functional::gradcorr(et2,vt2,v2,&chr,&rhopw,&ucell,stress2,false);
            XC_Functional::gradcorr(et2,vt2,v2,&chr,&rhopw,&ucell,stress2,true);

            GlobalV::NSPIN = 4;
            GlobalV::DOMAG = true;
            XC_Functional::gradcorr(et4,vt4,v4,&chr,&rhopw,&ucell,stress4,false); 
        }
};

TEST_F(XCTest_GRADCORR, set_xc_type)
{
    double et1_ref = -0.02083356309, vt1_ref = 0.1433626281;
    std::vector<double> v1_ref = {0,0.02403590412,0.01672229351,0.01340429824,0.01141731056};
    std::vector<double> stress1_ref = {-0.02536030461,0,0,-0.02536030461,-0.02536030461,0,-0.02536030461,-0.02536030461,-0.02536030461};

    EXPECT_NEAR(et1_ref,et1,1.0e-8);
    EXPECT_NEAR(vt1_ref,vt1,1.0e-8);
    for(int i=0;i<5;i++)
    {
        EXPECT_NEAR(v1(0,i),v1_ref[i],1.0e-8);
    }
    for(int i=0;i<9;i++)
    {
        EXPECT_NEAR(stress1[i],stress1_ref[i],1.0e-8);
    }

    double et2_ref = -0.02334069902, vt2_ref = 0.1585216001;
    std::vector<double> v2_ref1 = {0,0.01705561346,0.01079116099,0.008088425437,0.006519533763};
    std::vector<double> v2_ref2 = {0,0.09418744142,0.0767446959,0.06744349422,0.0613488043};
    std::vector<double> stress2_ref = {-0.0280735975,0,0,-0.0280735975,-0.0280735975,0,-0.0280735975,-0.0280735975,-0.0280735975};

    EXPECT_NEAR(et2_ref,et2,1.0e-8);
    EXPECT_NEAR(vt2_ref,vt2,1.0e-8);
    for(int i=0;i<5;i++)
    {
        EXPECT_NEAR(v2(0,i),v2_ref1[i],1.0e-8);
        EXPECT_NEAR(v2(1,i),v2_ref2[i],1.0e-8);
    }
    for(int i=0;i<9;i++)
    {
        EXPECT_NEAR(stress2[i],stress2_ref[i],1.0e-8);
    }

    double et4_ref = -0.1443518167, vt4_ref = 0.4761829579;
    std::vector<double> v4_ref1 = {0,0.03249745531,0.02610023454,0.02290998159,0.02087122756};
    std::vector<double> v4_ref2 = {0,0.003217727553,0.00258430831,0.002268426198,0.002066559469};
    std::vector<double> v4_ref3 = {0,0.03217727553,0.0258430831,0.02268426198,0.02066559469};
    std::vector<double> v4_ref4 = {0,0.003217727553,0.00258430831,0.002268426198,0.002066559469};
    EXPECT_NEAR(et4_ref,et4,1.0e-8);
    EXPECT_NEAR(vt4_ref,vt4,1.0e-8);
    for(int i=0;i<5;i++)
    {
        EXPECT_NEAR(v4(0,i),v4_ref1[i],1.0e-8);
        EXPECT_NEAR(v4(1,i),v4_ref2[i],1.0e-8);
        EXPECT_NEAR(v4(2,i),v4_ref3[i],1.0e-8);
        EXPECT_NEAR(v4(3,i),v4_ref4[i],1.0e-8);
    }
}

class XCTest_GRADWFC : public testing::Test
{
    protected:

        std::complex<double> * grad = nullptr;
        ~XCTest_GRADWFC()
        {
            delete[] grad;
        }

        void SetUp()
        {
            ModulePW::PW_Basis_K rhopw;
            rhopw.npwk = new int[1];
            rhopw.npwk[0] = 5;
            rhopw.nmaxgr = 5;
            rhopw.nrxx = 5;

            std::complex<double> rhog[5];
            for (int i=0;i<5;i++)
            {
                rhog[i] = double(i);
            }
            double tpiba = 1;

            grad = new std::complex<double>[15];

            XC_Functional::grad_wfc(rhog, 0, grad, &rhopw, tpiba);
        }
};

TEST_F(XCTest_GRADWFC, set_xc_type)
{

    for (int j=0;j<3;j++)
    {
        for (int i=0;i<5;i++)
        {
            EXPECT_NEAR(grad[i+j*5].real(),double(i*(j+1)),1e-8);
            EXPECT_NEAR(grad[i+j*5].imag(),0,1e-8);
        }
    }
}