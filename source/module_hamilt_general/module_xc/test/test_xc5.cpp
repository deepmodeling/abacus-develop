#include "../xc_functional.h"
#include "gtest/gtest.h"
#include "../exx_info.h"
#include "xc3_mock.h"
#include "module_base/matrix.h"
#include "../../../module_base/parallel_reduce.h"

/************************************************
*  unit test of functionals
***********************************************/

// For more information of the functions, check the comment of xc_functional.h
// Three functions are tested:
// v_xc, the unified interface of LDA and GGA functionals
// v_xc_libxc, called by v_xc, when we use functionals from LIBXC
// v_xc_meta, unified interface of mGGA functionals

template<>
void Parallel_Reduce::reduce_pool<double>(double& object)
{
#ifdef __MPI
	double swap = object;
	MPI_Allreduce(&swap , &object , 1, MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD);
#endif
    return;
}

class XCTest_VXC : public testing::Test
{
    protected:
    
        double et1 = 0, vt1 = 0;
        ModuleBase::matrix v1;

        double et2 = 0, vt2 = 0;
        ModuleBase::matrix v2;

        double et4 = 0, vt4 = 0;
        ModuleBase::matrix v4;

        void SetUp()
        {
            ModulePW::PW_Basis rhopw;
            UnitCell ucell;
            Charge chr;
            
            rhopw.nrxx = 5;
            rhopw.npw = 5;
            rhopw.nmaxgr = 5;
            rhopw.gcar = new ModuleBase::Vector3<double> [5];
            rhopw.nxyz = 1;

            ucell.tpiba = 1;
            ucell.magnet.lsign_ = true;
            ucell.cal_ux();
            ucell.omega = 1;

            chr.rhopw = &(rhopw);
            chr.rho = new double*[4];
            chr.rho[0] = new double[5];
            chr.rho[1] = new double[5];
            chr.rho[2] = new double[5];
            chr.rho[3] = new double[5];
            chr.rhog = new complex<double>*[2];
            chr.rhog[0] = new complex<double>[5];
            chr.rhog[1] = new complex<double>[5];

            chr.rho_core = new double[5];
            chr.rhog_core = new complex<double>[5];

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

            XC_Functional::set_xc_type("PBE");

            GlobalV::NSPIN = 1;
            std::tuple<double, double, ModuleBase::matrix> etxc_vtxc_v
                = XC_Functional::v_xc(rhopw.nrxx,&chr,&ucell);
            et1 = std::get<0>(etxc_vtxc_v);
            vt1 = std::get<1>(etxc_vtxc_v);
            v1  = std::get<2>(etxc_vtxc_v);

            GlobalV::NSPIN = 2;
            etxc_vtxc_v
                = XC_Functional::v_xc(rhopw.nrxx,&chr,&ucell);
            et2 = std::get<0>(etxc_vtxc_v);
            vt2 = std::get<1>(etxc_vtxc_v);
            v2  = std::get<2>(etxc_vtxc_v);

            GlobalV::NSPIN = 4;
            GlobalV::DOMAG = true;
            etxc_vtxc_v
                = XC_Functional::v_xc(rhopw.nrxx,&chr,&ucell);
            et4 = std::get<0>(etxc_vtxc_v);
            vt4 = std::get<1>(etxc_vtxc_v);
            v4  = std::get<2>(etxc_vtxc_v);
        }
};

TEST_F(XCTest_VXC, set_xc_type)
{

    EXPECT_NEAR(et1,-22.58755058,1.0e-8);
    EXPECT_NEAR(vt1,-29.58544157,1.0e-8);
    EXPECT_NEAR(v1(0,0),0,1.0e-8);
    EXPECT_NEAR(v1(0,1),-2.10436858,1.0e-8);
    EXPECT_NEAR(v1(0,2),-2.635713084,1.0e-8);
    EXPECT_NEAR(v1(0,3),-3.005351752,1.0e-8);
    EXPECT_NEAR(v1(0,4),-3.298397892,1.0e-8);

    EXPECT_NEAR(et2,-28.97838368,1.0e-8);
    EXPECT_NEAR(vt2,-38.15420234,1.0e-8);
    EXPECT_NEAR(v2(0,0),0,1.0e-8);
    EXPECT_NEAR(v2(0,1),-2.560885436,1.0e-8);
    EXPECT_NEAR(v2(0,2),-3.219339115,1.0e-8);
    EXPECT_NEAR(v2(0,3),-3.678772816,1.0e-8);
    EXPECT_NEAR(v2(0,4),-4.043604077,1.0e-8);
    EXPECT_NEAR(v2(1,0),0,1.0e-8);
    EXPECT_NEAR(v2(1,1),-1.394281236,1.0e-8);
    EXPECT_NEAR(v2(1,2),-1.739033356,1.0e-8);
    EXPECT_NEAR(v2(1,3),-1.97506482,1.0e-8);
    EXPECT_NEAR(v2(1,4),-2.160374198,1.0e-8);

    EXPECT_NEAR(et4,-27.40098253,1.0e-8);
    EXPECT_NEAR(vt4,-35.81948838,1.0e-8);
    EXPECT_NEAR(v4(0,0),0,1.0e-8);
    EXPECT_NEAR(v4(0,1),-1.559604078,1.0e-8);
    EXPECT_NEAR(v4(0,2),-1.920028447,1.0e-8);
    EXPECT_NEAR(v4(0,3),-2.168396069,1.0e-8);
    EXPECT_NEAR(v4(0,4),-2.36419592,1.0e-8);
    EXPECT_NEAR(v4(1,0),0,1.0e-8);
    EXPECT_NEAR(v4(1,1),-0.09308179605,1.0e-8);
    EXPECT_NEAR(v4(1,2),-0.123132664,1.0e-8);
    EXPECT_NEAR(v4(1,3),-0.144332804,1.0e-8);
    EXPECT_NEAR(v4(1,4),-0.16127282,1.0e-8);
    EXPECT_NEAR(v4(2,0),0,1.0e-8);
    EXPECT_NEAR(v4(2,1),-0.9308179605,1.0e-8);
    EXPECT_NEAR(v4(2,2),-1.23132664,1.0e-8);
    EXPECT_NEAR(v4(2,3),-1.44332804,1.0e-8);
    EXPECT_NEAR(v4(2,4),-1.6127282,1.0e-8);
    EXPECT_NEAR(v4(3,0),0,1.0e-8);
    EXPECT_NEAR(v4(3,1),-0.09308179605,1.0e-8);
    EXPECT_NEAR(v4(3,2),-0.123132664,1.0e-8);
    EXPECT_NEAR(v4(3,3),-0.144332804,1.0e-8);
    EXPECT_NEAR(v4(3,4),-0.16127282,1.0e-8);

}

class XCTest_VXC_Libxc : public testing::Test
{
    protected:
    
        double et1 = 0, vt1 = 0;
        ModuleBase::matrix v1;

        double et2 = 0, vt2 = 0;
        ModuleBase::matrix v2;

        double et4 = 0, vt4 = 0;
        ModuleBase::matrix v4;

        void SetUp()
        {
            ModulePW::PW_Basis rhopw;
            UnitCell ucell;
            Charge chr;

            rhopw.nrxx = 5;
            rhopw.npw = 5;
            rhopw.nmaxgr = 5;
            rhopw.gcar = new ModuleBase::Vector3<double> [5];
            rhopw.nxyz = 1;

            ucell.tpiba = 1;
            ucell.magnet.lsign_ = true;
            ucell.cal_ux();
            ucell.omega = 1;

            chr.rhopw = &(rhopw);
            chr.rho = new double*[4];
            chr.rho[0] = new double[5];
            chr.rho[1] = new double[5];
            chr.rho[2] = new double[5];
            chr.rho[3] = new double[5];
            chr.rhog = new complex<double>*[2];
            chr.rhog[0] = new complex<double>[5];
            chr.rhog[1] = new complex<double>[5];

            chr.rho_core = new double[5];
            chr.rhog_core = new complex<double>[5];

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

            XC_Functional::set_xc_type("GGA_X_PBE+GGA_C_PBE");

            GlobalV::NSPIN = 1;
            std::tuple<double, double, ModuleBase::matrix> etxc_vtxc_v
                = XC_Functional::v_xc(rhopw.nrxx,&chr,&ucell);
            et1 = std::get<0>(etxc_vtxc_v);
            vt1 = std::get<1>(etxc_vtxc_v);
            v1  = std::get<2>(etxc_vtxc_v);

            GlobalV::NSPIN = 2;
            etxc_vtxc_v
                = XC_Functional::v_xc(rhopw.nrxx,&chr,&ucell);
            et2 = std::get<0>(etxc_vtxc_v);
            vt2 = std::get<1>(etxc_vtxc_v);
            v2  = std::get<2>(etxc_vtxc_v);

            GlobalV::NSPIN = 4;
            GlobalV::DOMAG = true;
            etxc_vtxc_v
                = XC_Functional::v_xc(rhopw.nrxx,&chr,&ucell);
            et4 = std::get<0>(etxc_vtxc_v);
            vt4 = std::get<1>(etxc_vtxc_v);
            v4  = std::get<2>(etxc_vtxc_v);
        }
};

TEST_F(XCTest_VXC_Libxc, set_xc_type)
{

    EXPECT_NEAR(et1,-22.58754423,1.0e-8);
    EXPECT_NEAR(vt1,-29.58543393,1.0e-8);
    EXPECT_NEAR(v1(0,0),0,1.0e-8);
    EXPECT_NEAR(v1(0,1),-2.104367948,1.0e-8);
    EXPECT_NEAR(v1(0,2),-2.635712365,1.0e-8);
    EXPECT_NEAR(v1(0,3),-3.005350979,1.0e-8);
    EXPECT_NEAR(v1(0,4),-3.298397079,1.0e-8);

    EXPECT_NEAR(et2,-28.97838189,1.0e-8);
    EXPECT_NEAR(vt2,-38.1541987,1.0e-8);
    EXPECT_NEAR(v2(0,0),0,1.0e-8);
    EXPECT_NEAR(v2(0,1),-2.560885532,1.0e-8);
    EXPECT_NEAR(v2(0,2),-3.219339294,1.0e-8);
    EXPECT_NEAR(v2(0,3),-3.678773042,1.0e-8);
    EXPECT_NEAR(v2(0,4),-4.043604335,1.0e-8);
    EXPECT_NEAR(v2(1,0),0,1.0e-8);
    EXPECT_NEAR(v2(1,1),-1.394276473,1.0e-8);
    EXPECT_NEAR(v2(1,2),-1.739027899,1.0e-8);
    EXPECT_NEAR(v2(1,3),-1.975058937,1.0e-8);
    EXPECT_NEAR(v2(1,4),-2.160368003,1.0e-8);

    EXPECT_NEAR(et4,-27.28201062,1.0e-8);
    EXPECT_NEAR(vt4,-35.98253991,1.0e-8);
    EXPECT_NEAR(v4(0,0),0,1.0e-8);
    EXPECT_NEAR(v4(0,1),-1.268278149,1.0e-8);
    EXPECT_NEAR(v4(0,2),-1.598108222,1.0e-8);
    EXPECT_NEAR(v4(0,3),-1.828079634,1.0e-8);
    EXPECT_NEAR(v4(0,4),-2.010634115,1.0e-8);
    EXPECT_NEAR(v4(1,0),0,1.0e-8);
    EXPECT_NEAR(v4(1,1),-0.1255782493,1.0e-8);
    EXPECT_NEAR(v4(1,2),-0.1582362929,1.0e-8);
    EXPECT_NEAR(v4(1,3),-0.1810068558,1.0e-8);
    EXPECT_NEAR(v4(1,4),-0.1990824429,1.0e-8);
    EXPECT_NEAR(v4(2,0),0,1.0e-8);
    EXPECT_NEAR(v4(2,1),-1.255782493,1.0e-8);
    EXPECT_NEAR(v4(2,2),-1.582362929,1.0e-8);
    EXPECT_NEAR(v4(2,3),-1.810068558,1.0e-8);
    EXPECT_NEAR(v4(2,4),-1.990824429,1.0e-8);
    EXPECT_NEAR(v4(3,0),0,1.0e-8);
    EXPECT_NEAR(v4(3,1),-0.1255782493,1.0e-8);
    EXPECT_NEAR(v4(3,2),-0.1582362929,1.0e-8);
    EXPECT_NEAR(v4(3,3),-0.1810068558,1.0e-8);
    EXPECT_NEAR(v4(3,4),-0.1990824429,1.0e-8);
}

class XCTest_VXC_meta : public testing::Test
{
    protected:
    
        double et1 = 0, vt1 = 0;
        ModuleBase::matrix v1,vtau1;

        double et2 = 0, vt2 = 0;
        ModuleBase::matrix v2,vtau2;

        void SetUp()
        {
            ModulePW::PW_Basis rhopw;
            UnitCell ucell;
            Charge chr;

            rhopw.nrxx = 5;
            rhopw.npw = 5;
            rhopw.nmaxgr = 5;
            rhopw.gcar = new ModuleBase::Vector3<double> [5];
            rhopw.nxyz = 1;

            ucell.tpiba = 1;
            ucell.magnet.lsign_ = true;
            ucell.cal_ux();
            ucell.omega = 1;

            chr.rhopw = &(rhopw);
            chr.rho = new double*[2];
            chr.rho[0] = new double[5];
            chr.rho[1] = new double[5];
            chr.rhog = new complex<double>*[2];
            chr.rhog[0] = new complex<double>[5];
            chr.rhog[1] = new complex<double>[5];

            chr.rho_core = new double[5];
            chr.rhog_core = new complex<double>[5];

            for(int i=0;i<5;i++)
            {
                chr.rho[0][i] = double(i);
                chr.rho[1][i] = 0.1*double(i);
                chr.rhog[0][i] = chr.rho[0][i];
                chr.rhog[1][i] = chr.rho[1][i];
                chr.rho_core[i] = 0;
                chr.rhog_core[i] = 0;
                rhopw.gcar[i]= 1;
            }

            chr.kin_r = new double*[2];
            chr.kin_r[0] = new double[5];
            chr.kin_r[1] = new double[5];
            chr.kin_r[0][0] = 0;
            chr.kin_r[0][1] = 0.02403590412;
            chr.kin_r[0][2] = 0.01672229351;
            chr.kin_r[0][3] = 0.01340429824;
            chr.kin_r[0][4] = 0.01141731056;
            chr.kin_r[1][0] = 0.5;
            chr.kin_r[1][1] = 0.52403590412;
            chr.kin_r[1][2] = 0.51672229351;
            chr.kin_r[1][3] = 0.51340429824;
            chr.kin_r[1][4] = 0.51141731056;

            XC_Functional::set_xc_type("SCAN");

            GlobalV::NSPIN = 1;
            std::tuple<double, double, ModuleBase::matrix, ModuleBase::matrix> etxc_vtxc_v
                = XC_Functional::v_xc_meta(rhopw.nrxx,ucell.omega,ucell.tpiba,&chr);
            et1 = std::get<0>(etxc_vtxc_v);
            vt1 = std::get<1>(etxc_vtxc_v);
            v1  = std::get<2>(etxc_vtxc_v);
            vtau1 = std::get<3>(etxc_vtxc_v);

            GlobalV::NSPIN = 2;
            etxc_vtxc_v
                = XC_Functional::v_xc_meta(rhopw.nrxx,ucell.omega,ucell.tpiba,&chr);
            et2 = std::get<0>(etxc_vtxc_v);
            vt2 = std::get<1>(etxc_vtxc_v);
            v2  = std::get<2>(etxc_vtxc_v);
            vtau2 = std::get<3>(etxc_vtxc_v);
        }
};

TEST_F(XCTest_VXC_meta, set_xc_type)
{

    EXPECT_NEAR(et1,-25.13065363,1.0e-8);
    EXPECT_NEAR(vt1,-33.13880774,1.0e-8);
    EXPECT_NEAR(v1(0,0),0,1.0e-8);
    EXPECT_NEAR(v1(0,1),-2.336719556,1.0e-8);
    EXPECT_NEAR(v1(0,2),-2.942649664,1.0e-8);
    EXPECT_NEAR(v1(0,3),-3.36679035,1.0e-8);
    EXPECT_NEAR(v1(0,4),-3.704104452,1.0e-8);
    EXPECT_NEAR(vtau1(0,0),0,1.0e-8);
    EXPECT_NEAR(vtau1(0,1),0.0187099814,1.0e-8);
    EXPECT_NEAR(vtau1(0,2),0.01578002561,1.0e-8);
    EXPECT_NEAR(vtau1(0,3),0.01423896928,1.0e-8);
    EXPECT_NEAR(vtau1(0,4),0.01321861589,1.0e-8);

    EXPECT_NEAR(et2,-32.72218711,1.0e-8);
    EXPECT_NEAR(vt2,-43.31358017,1.0e-8);
    EXPECT_NEAR(v2(0,0),0,1.0e-8);
    EXPECT_NEAR(v2(0,1),-2.901190807,1.0e-8);
    EXPECT_NEAR(v2(0,2),-3.662642983,1.0e-8);
    EXPECT_NEAR(v2(0,3),-4.196098173,1.0e-8);
    EXPECT_NEAR(v2(0,4),-4.620454375,1.0e-8);
    EXPECT_NEAR(v2(1,0),0,1.0e-8);
    EXPECT_NEAR(v2(1,1),-1.285513329,1.0e-8);
    EXPECT_NEAR(v2(1,2),-1.795172177,1.0e-8);
    EXPECT_NEAR(v2(1,3),-2.064035864,1.0e-8);
    EXPECT_NEAR(v2(1,4),-2.275487119,1.0e-8);
    EXPECT_NEAR(vtau2(0,0),0,1.0e-8);
    EXPECT_NEAR(vtau2(0,1),0.01677946177,1.0e-8);
    EXPECT_NEAR(vtau2(0,2),0.01410816304,1.0e-8);
    EXPECT_NEAR(vtau2(0,3),0.01263339482,1.0e-8);
    EXPECT_NEAR(vtau2(0,4),0.01165715023,1.0e-8);
    EXPECT_NEAR(vtau2(1,0),0,1.0e-8);
    EXPECT_NEAR(vtau2(1,1),0.01591158497,1.0e-8);
    EXPECT_NEAR(vtau2(1,2),0.07990709956,1.0e-8);
    EXPECT_NEAR(vtau2(1,3),0.04145463825,1.0e-8);
    EXPECT_NEAR(vtau2(1,4),0.0311787189,1.0e-8);    
}


int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}