#include "../xc_functional.h"
#include "../libxc_abacus.h"
#include "gtest/gtest.h"
#include "xctest.h"
#include "../exx_info.h"
#include "xc3_mock.h"
#include "source_base/matrix.h"
#include "../../../source_base/parallel_reduce.h"

/************************************************
*  unit test of functionals
***********************************************/

// For more information of the functions, check the comment of xc_functional.h
// Three functions are tested:
// v_xc, the unified interface of LDA and GGA functionals
// v_xc_libxc, called by v_xc, when we use functionals from LIBXC
// v_xc_meta, unified interface of mGGA functionals

class XCTest_VXC : public XCTest
{
    protected:

        double et1 = 0, vt1 = 0;
        ModuleBase::matrix v1;

        double et2 = 0, vt2 = 0;
        ModuleBase::matrix v2;

        void SetUp()
        {
            // Define variables for parameters
            int nspin1 = 1;
            int nspin2 = 2;
            bool domag = false;
            bool domag_z = false;

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
            elecstate::cal_ux(ucell, 4);
            ucell.omega = 1;

            chr.rhopw = &(rhopw);
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

            XC_Functional::set_xc_type("PBE");

            const double hybrid_alpha = XC_Functional::get_hybrid_alpha();
            const double hse_omega = XC_Functional::get_hse_omega();
            std::tuple<double, double, ModuleBase::matrix> etxc_vtxc_v
                = XC_Functional::v_xc(rhopw.nrxx,&chr,&ucell,nspin1,domag,domag_z,0, hybrid_alpha, hse_omega);
            et1 = std::get<0>(etxc_vtxc_v);
            vt1 = std::get<1>(etxc_vtxc_v);
            v1  = std::get<2>(etxc_vtxc_v);

            etxc_vtxc_v
                = XC_Functional::v_xc(rhopw.nrxx,&chr,&ucell,nspin2,domag,domag_z,0, hybrid_alpha, hse_omega);
            et2 = std::get<0>(etxc_vtxc_v);
            vt2 = std::get<1>(etxc_vtxc_v);
            v2  = std::get<2>(etxc_vtxc_v);

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

}

class XCTest_VXC_Libxc : public XCTest
{
    protected:

        double et1 = 0, vt1 = 0;
        ModuleBase::matrix v1;

        double et2 = 0, vt2 = 0;
        ModuleBase::matrix v2;

        void SetUp()
        {
            // Define variables for parameters
            int nspin1 = 1;
            int nspin2 = 2;
            bool domag = false;
            bool domag_z = false;

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
            elecstate::cal_ux(ucell, 4);
            ucell.omega = 1;

            chr.rhopw = &(rhopw);
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

            XC_Functional::set_xc_type("GGA_X_PBE+GGA_C_PBE");

            const double hybrid_alpha = XC_Functional::get_hybrid_alpha();
            const double hse_omega = XC_Functional::get_hse_omega();
            std::tuple<double, double, ModuleBase::matrix> etxc_vtxc_v
                = XC_Functional::v_xc(rhopw.nrxx,&chr,&ucell,nspin1,domag,domag_z,0, hybrid_alpha, hse_omega);
            et1 = std::get<0>(etxc_vtxc_v);
            vt1 = std::get<1>(etxc_vtxc_v);
            v1  = std::get<2>(etxc_vtxc_v);

            etxc_vtxc_v
                = XC_Functional::v_xc(rhopw.nrxx,&chr,&ucell,nspin2,domag,domag_z,0, hybrid_alpha, hse_omega);
            et2 = std::get<0>(etxc_vtxc_v);
            vt2 = std::get<1>(etxc_vtxc_v);
            v2  = std::get<2>(etxc_vtxc_v);

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

}

class XCTest_VXC_meta : public XCTest
{
    protected:

        double et1 = 0, vt1 = 0;
        ModuleBase::matrix v1,vtau1;

        double et2 = 0, vt2 = 0;
        ModuleBase::matrix v2,vtau2;

        void SetUp()
        {
            // Define variables for parameters
            int nspin1 = 1;
            int nspin2 = 2;

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
            elecstate::cal_ux(ucell, 4);
            ucell.omega = 1;

            chr.rhopw = &(rhopw);
            chr.rho = new double*[2];
            chr.rho[0] = new double[5];
            chr.rho[1] = new double[5];
            chr.rhog = new std::complex<double>*[2];
            chr.rhog[0] = new std::complex<double>[5];
            chr.rhog[1] = new std::complex<double>[5];

            chr.rho_core = new double[5];
            chr.rhog_core = new std::complex<double>[5];

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

            const double hybrid_alpha = XC_Functional::get_hybrid_alpha();
            const double hse_omega = XC_Functional::get_hse_omega();
            std::tuple<double, double, ModuleBase::matrix, ModuleBase::matrix> etxc_vtxc_v
                = XC_Functional_Libxc::v_xc_meta(XC_Functional::get_func_id(), rhopw.nrxx,ucell.omega,ucell.tpiba,&chr,nspin1, hybrid_alpha, hse_omega);
            et1 = std::get<0>(etxc_vtxc_v);
            vt1 = std::get<1>(etxc_vtxc_v);
            v1  = std::get<2>(etxc_vtxc_v);
            vtau1 = std::get<3>(etxc_vtxc_v);

            etxc_vtxc_v
                = XC_Functional_Libxc::v_xc_meta(XC_Functional::get_func_id(), rhopw.nrxx,ucell.omega,ucell.tpiba,&chr,nspin2, hybrid_alpha, hse_omega);
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

/************************************************
 *  unit tests for the gga_grad keyword (nspin=4
 *  noncollinear GGA gradient methods)
 *
 *  Methods 2 and 3 share the same chain-rule spin-up/down gradients
 *    grad(rho_up/dn) = (grad(rho) +/- m_hat . grad(m))/2, m_hat = m/|m|
 *  but differ in the divergence of h = df/d(grad rho):
 *    gga_grad=2 (projected): v_mu -= m_hat_mu * div((h_up - h_dn)/2)
 *    gga_grad=3 (full SF):   v_mu -= div((h_up - h_dn)/2 * m_hat_mu)
 *  The two are identical when grad(m_hat) = 0 (magnetization direction
 *  uniform in space) and differ otherwise.
 ***********************************************/

namespace
{
constexpr int gga_grad_nrxx = 5;

// build a mock 4-component charge on the mocked 5-point grid.
// rotating=false: m = (0,0,mz), mz>0, so m_hat = (0,0,1) everywhere
// rotating=true:  m direction varies from point to point
struct Ns4Charge
{
    ModulePW::PW_Basis rhopw;
    UnitCell ucell;
    Charge chr;

    Ns4Charge(const bool rotating)
    {
        rhopw.nrxx = gga_grad_nrxx;
        rhopw.npw = gga_grad_nrxx;
        rhopw.nmaxgr = gga_grad_nrxx;
        rhopw.gcar = new ModuleBase::Vector3<double>[gga_grad_nrxx];
        rhopw.nxyz = 1;

        ucell.tpiba = 1;
        ucell.omega = 1;
        ucell.magnet.lsign_ = true;
        elecstate::cal_ux(ucell, 4);

        chr.rhopw = &(rhopw);
        chr.rho = new double*[4];
        for (int is = 0; is < 4; ++is)
        {
            chr.rho[is] = new double[gga_grad_nrxx];
        }
        chr.rhog = new std::complex<double>*[2];
        chr.rhog[0] = new std::complex<double>[gga_grad_nrxx];
        chr.rhog[1] = new std::complex<double>[gga_grad_nrxx];
        chr.rho_core = new double[gga_grad_nrxx];
        chr.rhog_core = new std::complex<double>[gga_grad_nrxx];

        for (int i = 0; i < gga_grad_nrxx; ++i)
        {
            chr.rho[0][i] = 2.0 + i;
            if (rotating)
            {
                chr.rho[1][i] = 0.10 * (i + 1);
                chr.rho[2][i] = 0.05 * (gga_grad_nrxx - i);
                chr.rho[3][i] = 0.20 * (i + 1);
            }
            else
            {
                chr.rho[1][i] = 0.0;
                chr.rho[2][i] = 0.0;
                chr.rho[3][i] = 0.2 * (i + 1);
            }
            chr.rhog[0][i] = chr.rho[0][i];
            chr.rhog[1][i] = chr.rho[1][i];
            chr.rho_core[i] = 0;
            chr.rhog_core[i] = 0;
            rhopw.gcar[i] = 1;
        }
    }
};

// run XC_Functional::v_xc for nspin=4 with noncollinear magnetism
std::tuple<double, double, ModuleBase::matrix> run_vxc_nspin4(
    const std::string& functional,
    const bool rotating,
    const int gga_grad)
{
    Ns4Charge mock(rotating);
    XC_Functional::set_xc_type(functional);
    return XC_Functional::v_xc(gga_grad_nrxx,
                               &mock.chr,
                               &mock.ucell,
                               4,
                               true,
                               false,
                               gga_grad,
                               XC_Functional::get_hybrid_alpha(),
                               XC_Functional::get_hse_omega());
}

// compare two (etxc, vtxc, v) results
void expect_vxc_equal(const std::tuple<double, double, ModuleBase::matrix>& a,
                      const std::tuple<double, double, ModuleBase::matrix>& b,
                      const double tol)
{
    EXPECT_NEAR(std::get<0>(a), std::get<0>(b), tol);
    EXPECT_NEAR(std::get<1>(a), std::get<1>(b), tol);
    const ModuleBase::matrix& va = std::get<2>(a);
    const ModuleBase::matrix& vb = std::get<2>(b);
    ASSERT_EQ(va.nr, vb.nr);
    ASSERT_EQ(va.nc, vb.nc);
    for (int ir = 0; ir < va.nr; ++ir)
    {
        for (int ic = 0; ic < va.nc; ++ic)
        {
            EXPECT_NEAR(va(ir, ic), vb(ir, ic), tol);
        }
    }
}
} // namespace

// m_hat = m/|m|, zero where |m| ~ 0
TEST(GgaGradTools, ComputeMagPartNspin4)
{
    Charge chr;
    chr.rho = new double*[4];
    for (int is = 0; is < 4; ++is)
    {
        chr.rho[is] = new double[3];
    }
    chr.rho[1][0] = 3.0; chr.rho[2][0] = 0.0; chr.rho[3][0] = 4.0; // |m|=5
    chr.rho[1][1] = 0.0; chr.rho[2][1] = 0.0; chr.rho[3][1] = 0.0; // |m|=0
    chr.rho[1][2] = 1.0; chr.rho[2][2] = 2.0; chr.rho[3][2] = 2.0; // |m|=3

    const std::vector<double> mp = XC_Functional_Libxc::compute_mag_part_nspin4(3, &chr);

    EXPECT_NEAR(mp[0], 3.0 / 5.0, 1e-15);
    EXPECT_NEAR(mp[3], 0.0, 1e-15);
    EXPECT_NEAR(mp[6], 4.0 / 5.0, 1e-15);
    EXPECT_NEAR(mp[1], 0.0, 1e-15);
    EXPECT_NEAR(mp[4], 0.0, 1e-15);
    EXPECT_NEAR(mp[7], 0.0, 1e-15);
    EXPECT_NEAR(mp[2], 1.0 / 3.0, 1e-15);
    EXPECT_NEAR(mp[5], 2.0 / 3.0, 1e-15);
    EXPECT_NEAR(mp[8], 2.0 / 3.0, 1e-15);
}

// v_tot = 0.5*(v_up+v_dn), v_mu = 0.5*(v_up-v_dn)*m_hat_mu
TEST(GgaGradTools, ConvertVNspin4Sf)
{
    Ns4Charge mock(false); // m_hat = (0,0,1)
    const std::vector<double> mag_part
        = XC_Functional_Libxc::compute_mag_part_nspin4(gga_grad_nrxx, &mock.chr);

    ModuleBase::matrix v(2, gga_grad_nrxx);
    for (int ir = 0; ir < gga_grad_nrxx; ++ir)
    {
        v(0, ir) = 1.0 + ir;
        v(1, ir) = 0.5 * ir;
    }
    const ModuleBase::matrix v4
        = XC_Functional_Libxc::convert_v_nspin4_sf(gga_grad_nrxx, &mock.chr, mag_part, v);

    for (int ir = 0; ir < gga_grad_nrxx; ++ir)
    {
        EXPECT_NEAR(v4(0, ir), 0.5 * (v(0, ir) + v(1, ir)), 1e-14);
        EXPECT_NEAR(v4(1, ir), 0.0, 1e-14);
        EXPECT_NEAR(v4(2, ir), 0.0, 1e-14);
        EXPECT_NEAR(v4(3, ir), 0.5 * (v(0, ir) - v(1, ir)), 1e-14);
    }
}

// original conversion: has_mag=false leaves magnetic channels zero
TEST(GgaGradTools, ConvertVNspin4HasMag)
{
    Ns4Charge mock(false);
    std::vector<double> amag(gga_grad_nrxx);
    for (int ir = 0; ir < gga_grad_nrxx; ++ir)
    {
        amag[ir] = mock.chr.rho[3][ir];
    }
    ModuleBase::matrix v(2, gga_grad_nrxx);
    for (int ir = 0; ir < gga_grad_nrxx; ++ir)
    {
        v(0, ir) = 1.0 + ir;
        v(1, ir) = 0.5 * ir;
    }

    const ModuleBase::matrix v_nomag
        = XC_Functional_Libxc::convert_v_nspin4(gga_grad_nrxx, &mock.chr, amag, v, false);
    for (int ir = 0; ir < gga_grad_nrxx; ++ir)
    {
        EXPECT_NEAR(v_nomag(0, ir), 0.5 * (v(0, ir) + v(1, ir)), 1e-14);
        EXPECT_NEAR(v_nomag(1, ir), 0.0, 1e-14);
        EXPECT_NEAR(v_nomag(2, ir), 0.0, 1e-14);
        EXPECT_NEAR(v_nomag(3, ir), 0.0, 1e-14);
    }

    const ModuleBase::matrix v_mag
        = XC_Functional_Libxc::convert_v_nspin4(gga_grad_nrxx, &mock.chr, amag, v, true);
    for (int ir = 0; ir < gga_grad_nrxx; ++ir)
    {
        const double vs = 0.5 * (v(0, ir) - v(1, ir));
        EXPECT_NEAR(v_mag(3, ir), vs * mock.chr.rho[3][ir] / amag[ir], 1e-14);
    }
}

class GgaGradDh : public testing::Test
{
  protected:
    // dh from cal_dh_sf for gga_grad=2 and 3 on a given magnetization pattern
    void run(const bool rotating,
             std::vector<std::vector<double>>& dh2,
             std::vector<std::vector<double>>& dh3,
             std::vector<double>& mag_part)
    {
        Ns4Charge mock(rotating);
        mag_part = XC_Functional_Libxc::compute_mag_part_nspin4(gga_grad_nrxx, &mock.chr);

        const std::tuple<std::vector<double>, std::vector<double>> rho_amag
            = XC_Functional_Libxc::convert_rho_amag_nspin4(2, gga_grad_nrxx, &mock.chr);
        const std::vector<double>& rho = std::get<0>(rho_amag);
        const std::vector<std::vector<ModuleBase::Vector3<double>>> gdr
            = XC_Functional_Libxc::cal_gdr_sf(2, gga_grad_nrxx, rho, mag_part, mock.ucell.tpiba, &mock.chr);

        std::vector<double> sgn(gga_grad_nrxx * 2, 1.0);
        std::vector<double> vsigma(gga_grad_nrxx * 3);
        for (int ir = 0; ir < gga_grad_nrxx; ++ir)
        {
            for (int j = 0; j < 3; ++j)
            {
                vsigma[ir * 3 + j] = 0.2 + 0.1 * ir + 0.05 * j;
            }
        }

        dh2 = XC_Functional_Libxc::cal_dh_sf(
            2, gga_grad_nrxx, sgn, gdr, vsigma, mag_part, 2, mock.ucell.tpiba, &mock.chr);
        dh3 = XC_Functional_Libxc::cal_dh_sf(
            2, gga_grad_nrxx, sgn, gdr, vsigma, mag_part, 3, mock.ucell.tpiba, &mock.chr);
    }
};

// uniform m_hat => grad(m_hat) = 0 => projected and full SF divergences agree
TEST_F(GgaGradDh, UniformDirectionMethodsAgree)
{
    std::vector<std::vector<double>> dh2, dh3;
    std::vector<double> mag_part;
    run(false, dh2, dh3, mag_part);

    for (int is = 0; is < 4; ++is)
    {
        for (int ir = 0; ir < gga_grad_nrxx; ++ir)
        {
            EXPECT_NEAR(dh2[is][ir], dh3[is][ir], 1e-12);
        }
    }
}

// for gga_grad=2, dh_mu = m_hat_mu * div((h_up-h_dn)/2), so the magnetic
// channels satisfy dh_mu = m_hat_mu * (m_hat . dh) exactly
TEST_F(GgaGradDh, ProjectedDivergenceIsProjection)
{
    std::vector<std::vector<double>> dh2, dh3;
    std::vector<double> mag_part;
    run(true, dh2, dh3, mag_part);

    for (int ir = 0; ir < gga_grad_nrxx; ++ir)
    {
        double proj = 0.0;
        for (int mu = 1; mu < 4; ++mu)
        {
            proj += mag_part[ir + (mu - 1) * gga_grad_nrxx] * dh2[mu][ir];
        }
        for (int mu = 1; mu < 4; ++mu)
        {
            EXPECT_NEAR(dh2[mu][ir], mag_part[ir + (mu - 1) * gga_grad_nrxx] * proj, 1e-12);
        }
    }
}

// NOTE on the mocked FFT: the mock derivative is pointwise
// (grad(f)[ir] ~ f[ir], div(h)[ir] ~ h[ir]), so multiplication by a
// scalar field commutes with the divergence and the projected (2) and
// full (3) SF divergences coincide under this mock. The two methods can
// only be distinguished with a real (nonlocal) FFT, e.g. in integration
// tests. The tests below therefore anchor the wiring (exact values,
// projection identities, crash-free dispatch) rather than the 2-vs-3
// numerical difference.

// built-in functionals: v_xc dispatches nspin=4 + gga_grad=2/3 to the SF builtin
TEST(GgaGradVxc, BuiltinUniformDirectionMethodsAgree)
{
    const auto r2 = run_vxc_nspin4("PBE", false, 2);
    const auto r3 = run_vxc_nspin4("PBE", false, 3);
    EXPECT_EQ(std::get<2>(r2).nr, 4);
    expect_vxc_equal(r2, r3, 1e-10);
}

// gga_grad=0 keeps the original built-in algorithm and must not crash
TEST(GgaGradVxc, BuiltinOriginalAlgorithmRuns)
{
    const auto r0 = run_vxc_nspin4("PBE", true, 0);
    EXPECT_EQ(std::get<2>(r0).nr, 4);
    EXPECT_TRUE(std::isfinite(std::get<0>(r0)));
    EXPECT_TRUE(std::isfinite(std::get<1>(r0)));
}

// regression anchors for the built-in SF path (gga_grad=3)
TEST(GgaGradVxc, BuiltinSfAnchoredValues)
{
    const auto r = run_vxc_nspin4("PBE", true, 3);
    const ModuleBase::matrix& v = std::get<2>(r);
    EXPECT_NEAR(std::get<0>(r), -5.1906253324e+01, 1.0e-8);
    EXPECT_NEAR(std::get<1>(r), -6.8184946486e+01, 1.0e-8);
    EXPECT_NEAR(v(0, 0), -2.6284212276e+00, 1.0e-8);
    EXPECT_NEAR(v(0, 4), -3.7476729308e+00, 1.0e-8);
    EXPECT_NEAR(v(1, 0), -3.7774725540e-02, 1.0e-8);
    EXPECT_NEAR(v(1, 4), -9.2204398939e-02, 1.0e-8);
    EXPECT_NEAR(v(2, 0), -9.4436813851e-02, 1.0e-8);
    EXPECT_NEAR(v(2, 4), -9.2204398939e-03, 1.0e-8);
    EXPECT_NEAR(v(3, 0), -7.5549451081e-02, 1.0e-8);
    EXPECT_NEAR(v(3, 4), -1.8440879788e-01, 1.0e-8);
}

// LIBXC functionals: gga_grad=2/3 select the SF path in v_xc_libxc
TEST(GgaGradVxc, LibxcUniformDirectionMethodsAgree)
{
    const auto r2 = run_vxc_nspin4("GGA_X_PBE+GGA_C_PBE", false, 2);
    const auto r3 = run_vxc_nspin4("GGA_X_PBE+GGA_C_PBE", false, 3);
    EXPECT_EQ(std::get<2>(r2).nr, 4);
    expect_vxc_equal(r2, r3, 1e-10);
}

// for LIBXC, gga_grad=0 and 1 both keep the original collinear algorithm
TEST(GgaGradVxc, LibxcZeroEqualsOne)
{
    const auto r0 = run_vxc_nspin4("GGA_X_PBE+GGA_C_PBE", true, 0);
    const auto r1 = run_vxc_nspin4("GGA_X_PBE+GGA_C_PBE", true, 1);
    expect_vxc_equal(r0, r1, 1e-12);
}

// regression anchors for the LIBXC SF path (gga_grad=3); this path goes
// through the SF branch of convert_vtxc_v
TEST(GgaGradVxc, LibxcSfAnchoredValues)
{
    const auto r = run_vxc_nspin4("GGA_X_PBE+GGA_C_PBE", true, 3);
    const ModuleBase::matrix& v = std::get<2>(r);
    EXPECT_NEAR(std::get<0>(r), -5.1906239921e+01, 1.0e-8);
    EXPECT_NEAR(std::get<1>(r), -6.8184928281e+01, 1.0e-8);
    EXPECT_NEAR(v(0, 0), -2.6284203722e+00, 1.0e-8);
    EXPECT_NEAR(v(0, 4), -3.7476719779e+00, 1.0e-8);
    EXPECT_NEAR(v(1, 0), -3.7774739650e-02, 1.0e-8);
    EXPECT_NEAR(v(1, 4), -9.2204427285e-02, 1.0e-8);
    EXPECT_NEAR(v(2, 0), -9.4436849124e-02, 1.0e-8);
    EXPECT_NEAR(v(2, 4), -9.2204427285e-03, 1.0e-8);
    EXPECT_NEAR(v(3, 0), -7.5549479300e-02, 1.0e-8);
    EXPECT_NEAR(v(3, 4), -1.8440885457e-01, 1.0e-8);
}


int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}