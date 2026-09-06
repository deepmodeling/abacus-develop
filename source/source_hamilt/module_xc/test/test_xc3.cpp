#include "gtest/gtest.h"
#include "xctest.h"
#include "../xc_functional.h"
#include "../exx_info.h"
#include "xc3_mock.h"
#include "source_base/matrix.h"
#include "source_cell/cal_ux.h"

/************************************************
*  unit test of functionals
***********************************************/

// For more information of the functions, check the comment of xc_functional.h
// Two functions are tested:
// gradcorr, which calculates the gradient part of GGA functional
// gradwfc, which is used to obtain the derivative of wavefunction



class XCTest_GRADCORR : public XCTest
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
            // Define variables for parameters
            int nspin1 = 1;
            int nspin2 = 2;
            int nspin4 = 4;
            bool domag = false;
            bool domag_z = false;
            bool domag_true = true;

            ModulePW::PW_Basis rhopw;
            UnitCell ucell;
            Charge chr;

            rhopw.nrxx = 5;
            rhopw.npw = 5;
            rhopw.nmaxgr = 5;
            rhopw.gcar = new ModuleBase::Vector3<double> [5];

            ucell.tpiba = 1;
            ucell.magnet.lsign_ = true;
            unitcell::cal_ux(ucell, 4);

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

            double hybrid_alpha = 0.0;
            double hse_omega = 0.0;
            XC_Functional::gradcorr(et1,vt1,v1,&chr,&rhopw,&ucell,stress1,false,nspin1,domag,domag_z, hybrid_alpha, hse_omega);
            XC_Functional::gradcorr(et1,vt1,v1,&chr,&rhopw,&ucell,stress1,true,nspin1,domag,domag_z, hybrid_alpha, hse_omega);

            XC_Functional::gradcorr(et2,vt2,v2,&chr,&rhopw,&ucell,stress2,false,nspin2,domag,domag_z, hybrid_alpha, hse_omega);
            XC_Functional::gradcorr(et2,vt2,v2,&chr,&rhopw,&ucell,stress2,true,nspin2,domag,domag_z, hybrid_alpha, hse_omega);

            XC_Functional::gradcorr(et4,vt4,v4,&chr,&rhopw,&ucell,stress4,false,nspin4,domag_true,domag_z, hybrid_alpha, hse_omega);
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

// Regression test for an out-of-bounds read in gradcorr.
// `set_xc_type("HF")` sets func_type = 4 (hybrid) but leaves func_id EMPTY, because pure
// Hartree-Fock has no (semi-)local functional. gradcorr used to index func_id[0]/func_id[1]
// unconditionally (see also the commented-out loop in XC_Functional::gcx_spin), reading back
// whatever the previously selected functional had left in the cleared-but-not-deallocated
// buffer. In the two-level EXX loop that predecessor is PBE (set_xc_first_loop), so a pure-HF
// run silently picked up a PBE gradient correction (issue deepmodeling/abacus-develop#5404)
//
// Setting PBE *before* HF is therefore an essential part of this test: it is what leaves the
// stale ids in the buffer, and without it the regression does not reproduce.
class XCTest_GRADCORR_HF : public XCTest
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
            const int nspin1 = 1;
            const int nspin2 = 2;
            const int nspin4 = 4;
            const bool domag = false;
            const bool domag_z = false;
            const bool domag_true = true;

            ModulePW::PW_Basis rhopw;
            UnitCell ucell;
            Charge chr;

            rhopw.nrxx = 5;
            rhopw.npw = 5;
            rhopw.nmaxgr = 5;
            rhopw.gcar = new ModuleBase::Vector3<double> [5];

            ucell.tpiba = 1;
            ucell.magnet.lsign_ = true;
            unitcell::cal_ux(ucell, 4);

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

            // leave PBE ids in func_id's buffer, then switch to pure HF
            XC_Functional::set_xc_type("PBE");
            XC_Functional::set_xc_type("HF");

            const double hybrid_alpha = 1.0;
            const double hse_omega = 0.0;
            XC_Functional::gradcorr(et1,vt1,v1,&chr,&rhopw,&ucell,stress1,false,nspin1,domag,domag_z, hybrid_alpha, hse_omega);
            XC_Functional::gradcorr(et1,vt1,v1,&chr,&rhopw,&ucell,stress1,true, nspin1,domag,domag_z, hybrid_alpha, hse_omega);

            XC_Functional::gradcorr(et2,vt2,v2,&chr,&rhopw,&ucell,stress2,false,nspin2,domag,domag_z, hybrid_alpha, hse_omega);
            XC_Functional::gradcorr(et2,vt2,v2,&chr,&rhopw,&ucell,stress2,true, nspin2,domag,domag_z, hybrid_alpha, hse_omega);

            XC_Functional::gradcorr(et4,vt4,v4,&chr,&rhopw,&ucell,stress4,false,nspin4,domag_true,domag_z, hybrid_alpha, hse_omega);
            XC_Functional::gradcorr(et4,vt4,v4,&chr,&rhopw,&ucell,stress4,true, nspin4,domag_true,domag_z, hybrid_alpha, hse_omega);
        }
};

TEST_F(XCTest_GRADCORR_HF, no_local_functional)
{
    // the state the fixture ran under: hybrid, but with no (semi-)local functional
    EXPECT_EQ(XC_Functional::get_func_type(), 4);
    EXPECT_TRUE(XC_Functional::get_func_id().empty());

    // pure HF must contribute exactly nothing through the gradient correction
    EXPECT_DOUBLE_EQ(et1, 0.0);
    EXPECT_DOUBLE_EQ(vt1, 0.0);
    EXPECT_DOUBLE_EQ(et2, 0.0);
    EXPECT_DOUBLE_EQ(vt2, 0.0);
    EXPECT_DOUBLE_EQ(et4, 0.0);
    EXPECT_DOUBLE_EQ(vt4, 0.0);

    for(int i=0;i<5;i++)
    {
        EXPECT_DOUBLE_EQ(v1(0,i), 0.0);
        EXPECT_DOUBLE_EQ(v2(0,i), 0.0);
        EXPECT_DOUBLE_EQ(v2(1,i), 0.0);
        for(int is=0;is<4;is++)
        {
            EXPECT_DOUBLE_EQ(v4(is,i), 0.0);
        }
    }

    // `Stress_Func::stress_gga` guards only on func_type, so it reaches this path for HF and
    // then reads stress_gga[0..8]: the tensor must come back sized and zeroed, not empty.
    ASSERT_EQ(stress1.size(), 9);
    ASSERT_EQ(stress2.size(), 9);
    ASSERT_EQ(stress4.size(), 9);
    for(int i=0;i<9;i++)
    {
        EXPECT_DOUBLE_EQ(stress1[i], 0.0);
        EXPECT_DOUBLE_EQ(stress2[i], 0.0);
        EXPECT_DOUBLE_EQ(stress4[i], 0.0);
    }
}

class XCTest_GRADWFC : public XCTest
{
    protected:

        std::complex<double> * grad = nullptr;
        ModuleBase::Vector3<double> * gcar_wrapper = nullptr;
        ModuleBase::Vector3<double> * kvec_c_wrapper = nullptr;
        ~XCTest_GRADWFC()
        {
            delete [] grad;
            delete [] gcar_wrapper;
            delete [] kvec_c_wrapper;
        }

        void SetUp()
        {
            ModulePW::PW_Basis_K rhopw;
            rhopw.npwk = new int[1];
            rhopw.npwk[0] = 5;
            rhopw.npwk_max = 5;
            rhopw.nmaxgr = 5;
            rhopw.nrxx = 5;
            rhopw.nks = 1;
            gcar_wrapper = new ModuleBase::Vector3<double>[rhopw.npwk_max];
            for (int ii = 0; ii < rhopw.npwk_max; ii++) {
                gcar_wrapper[ii] = ModuleBase::Vector3<double>(0,0,0);
}
            kvec_c_wrapper = new ModuleBase::Vector3<double>(1,2,3);

            rhopw.gcar = gcar_wrapper;
            rhopw.kvec_c = kvec_c_wrapper;

            std::complex<double> rhog[5];
            for (int i=0;i<5;i++)
            {
                rhog[i] = double(i);
            }
            double tpiba = 1;

            grad = new std::complex<double>[15];

            XC_Functional::grad_wfc<std::complex<double>, base_device::DEVICE_CPU>(0, tpiba, &rhopw, rhog, grad);
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