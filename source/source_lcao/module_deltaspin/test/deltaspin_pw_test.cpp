#include "gtest/gtest.h"
#include <complex>
#include <cmath>
#include <vector>

#define private public
#include "source_io/module_parameter/parameter.h"
#undef private

/***********************************************************************
 * Unit tests for DeltaSpin PW support
 *
 * Strategy: test the core arithmetic of calculate_delta_hcc and
 * cal_Mi_pw as pure formulas — no OnsiteProjector or full ABACUS
 * framework needed.
 ***********************************************************************/

class DeltaSpinPwTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}
};

// =====================================================================
// calculate_delta_hcc: ps array construction (npol=2, Pauli matrix)
// =====================================================================

TEST_F(DeltaSpinPwTest, DeltaHcc_Npol2_SingleAtom)
{
    // npol=2: for each (ib, ip):
    //   ps[becpind]      += coeff0 * becp1 + coeff2 * becp2
    //   ps[becpind+nkb]  += coeff1 * becp1 + coeff3 * becp2
    // where coeff0 = (lambda_z, 0), coeff1 = (lambda_x, lambda_y),
    //       coeff2 = (lambda_x, -lambda_y), coeff3 = (-lambda_z, 0)

    const int nat = 1;
    const int nproj = 2; // 2 projectors for this atom
    const int nbands = 1;
    const int nkb = nproj; // total projectors = nproj for single atom
    const int npol = 2;

    // delta_lambda for atom 0
    struct { double x, y, z; } delta_lambda = {0.5, 0.3, 0.8};

    const std::complex<double> coeff0(delta_lambda.z, 0.0);           // (0.8, 0)
    const std::complex<double> coeff1(delta_lambda.x, delta_lambda.y); // (0.5, 0.3)
    const std::complex<double> coeff2(delta_lambda.x, -delta_lambda.y);// (0.5, -0.3)
    const std::complex<double> coeff3(-delta_lambda.z, 0.0);          // (-0.8, 0)

    // becp: layout [ib * npol * nkb + sum + ip] for up, +nkb for down
    std::vector<std::complex<double>> becp(nbands * npol * nkb, {0.0, 0.0});
    // band 0, projector 0
    becp[0 * npol * nkb + 0] = {1.0, 0.2};       // becp_up[0]
    becp[0 * npol * nkb + 0 + nkb] = {0.3, -0.1}; // becp_dn[0]
    // band 0, projector 1
    becp[0 * npol * nkb + 1] = {0.5, 0.0};        // becp_up[1]
    becp[0 * npol * nkb + 1 + nkb] = {0.0, 0.7};  // becp_dn[1]

    std::vector<std::complex<double>> ps(nbands * npol * nkb, {0.0, 0.0});

    int sum = 0;
    for(int ib = 0; ib < nbands * npol; ib += npol)
    {
        for(int ip = 0; ip < nproj; ip++)
        {
            const int becpind = ib * nkb + sum + ip;
            const std::complex<double> becp1 = becp[becpind];
            const std::complex<double> becp2 = becp[becpind + nkb];
            ps[becpind] += coeff0 * becp1 + coeff2 * becp2;
            ps[becpind + nkb] += coeff1 * becp1 + coeff3 * becp2;
        }
    }

    // Verify projector 0:
    // ps_up[0] = (0.8,0)*(1.0,0.2) + (0.5,-0.3)*(0.3,-0.1)
    //          = (0.8, 0.16) + (0.15-0.03, -0.05-0.09) = (0.8,0.16) + (0.12,-0.14)
    //          = (0.92, 0.02)
    EXPECT_NEAR(ps[0].real(), 0.92, 1e-12);
    EXPECT_NEAR(ps[0].imag(), 0.02, 1e-12);

    // ps_dn[0] = (0.5,0.3)*(1.0,0.2) + (-0.8,0)*(0.3,-0.1)
    //          = (0.5-0.06, 0.3+0.1) + (-0.24, 0.08)
    //          = (0.44, 0.4) + (-0.24, 0.08) = (0.20, 0.48)
    EXPECT_NEAR(ps[0 + nkb].real(), 0.20, 1e-12);
    EXPECT_NEAR(ps[0 + nkb].imag(), 0.48, 1e-12);
}

TEST_F(DeltaSpinPwTest, DeltaHcc_Npol2_MultiAtom)
{
    // Two atoms: verify sum offset advances correctly
    const int nat = 2;
    const int nproj_0 = 1, nproj_1 = 1;
    const int nkb = nproj_0 + nproj_1; // 2
    const int nbands = 1;
    const int npol = 2;

    struct Vec3 { double x, y, z; };
    Vec3 delta_lambda[2] = {{1.0, 0.0, 0.0}, {0.0, 0.0, 2.0}};

    std::vector<std::complex<double>> becp(nbands * npol * nkb, {0.0, 0.0});
    // atom 0, proj 0: becp_up = (1,0), becp_dn = (0,0)
    becp[0] = {1.0, 0.0};
    becp[0 + nkb] = {0.0, 0.0};
    // atom 1, proj 0: becp_up = (0,0), becp_dn = (1,0)
    becp[1] = {0.0, 0.0};
    becp[1 + nkb] = {1.0, 0.0};

    std::vector<std::complex<double>> ps(nbands * npol * nkb, {0.0, 0.0});
    int nh_iat[2] = {nproj_0, nproj_1};

    int sum = 0;
    for(int iat = 0; iat < nat; iat++)
    {
        const std::complex<double> c0(delta_lambda[iat].z, 0.0);
        const std::complex<double> c1(delta_lambda[iat].x, delta_lambda[iat].y);
        const std::complex<double> c2(delta_lambda[iat].x, -delta_lambda[iat].y);
        const std::complex<double> c3(-delta_lambda[iat].z, 0.0);
        for(int ib = 0; ib < nbands * npol; ib += npol)
        {
            for(int ip = 0; ip < nh_iat[iat]; ip++)
            {
                const int becpind = ib * nkb + sum + ip;
                const std::complex<double> b1 = becp[becpind];
                const std::complex<double> b2 = becp[becpind + nkb];
                ps[becpind] += c0 * b1 + c2 * b2;
                ps[becpind + nkb] += c1 * b1 + c3 * b2;
            }
        }
        sum += nh_iat[iat];
    }

    // atom 0: lambda=(1,0,0), becp_up=(1,0), becp_dn=(0,0)
    // ps_up[0] = (0,0)*(1,0) + (1,0)*(0,0) = 0
    // ps_dn[0] = (1,0)*(1,0) + (0,0)*(0,0) = (1,0)
    EXPECT_NEAR(ps[0].real(), 0.0, 1e-12);
    EXPECT_NEAR(ps[0 + nkb].real(), 1.0, 1e-12);

    // atom 1: lambda=(0,0,2), becp_up=(0,0), becp_dn=(1,0)
    // ps_up[1] = (2,0)*(0,0) + (0,0)*(1,0) = 0
    // ps_dn[1] = (0,0)*(0,0) + (-2,0)*(1,0) = (-2,0)
    EXPECT_NEAR(ps[1].real(), 0.0, 1e-12);
    EXPECT_NEAR(ps[1 + nkb].real(), -2.0, 1e-12);
}

TEST_F(DeltaSpinPwTest, DeltaHcc_Npol1_SignPositive)
{
    // npol=1: ps[becpind] += sign * lambda_z * becp1
    const int nat = 1;
    const int nproj = 2;
    const int nkb = nproj;
    const int nbands = 1;
    const int sign = 1;
    const double lambda_z = 0.5;

    std::vector<std::complex<double>> becp(nbands * nkb, {0.0, 0.0});
    becp[0] = {1.0, 0.3};
    becp[1] = {0.0, -0.5};

    std::vector<std::complex<double>> ps(nbands * nkb, {0.0, 0.0});
    double coeff = lambda_z * sign;
    int sum = 0;
    for(int ib = 0; ib < nbands; ib++)
    {
        for(int ip = 0; ip < nproj; ip++)
        {
            const int becpind = ib * nkb + sum + ip;
            ps[becpind] += coeff * becp[becpind];
        }
    }

    // ps[0] = 0.5 * (1.0, 0.3) = (0.5, 0.15)
    EXPECT_NEAR(ps[0].real(), 0.5, 1e-12);
    EXPECT_NEAR(ps[0].imag(), 0.15, 1e-12);
    // ps[1] = 0.5 * (0, -0.5) = (0, -0.25)
    EXPECT_NEAR(ps[1].real(), 0.0, 1e-12);
    EXPECT_NEAR(ps[1].imag(), -0.25, 1e-12);
}

// =====================================================================
// cal_Mi_pw: magnetization accumulation from becp
// =====================================================================

TEST_F(DeltaSpinPwTest, MiPw_Npol1_SpinUp)
{
    // npol=1, nspin=2: Mi.z += sign * weight * |becp|^2
    // spin-up (sign=+1)
    const int nkb = 3;
    const int nbands = 2;
    const int sign = 1;
    const double weights[2] = {1.0, 0.5};

    std::vector<std::complex<double>> becp(nbands * nkb, {0.0, 0.0});
    // band 0
    becp[0] = {0.8, 0.0};
    becp[1] = {0.0, 0.6};
    becp[2] = {0.3, 0.4};
    // band 1
    becp[3] = {0.5, 0.0};
    becp[4] = {0.0, 0.0};
    becp[5] = {1.0, 0.0};

    // Single atom with nh=3
    double Mi_z = 0.0;
    for(int ib = 0; ib < nbands; ib++)
    {
        const double weight = weights[ib];
        double occ = 0.0;
        for(int ih = 0; ih < nkb; ih++)
        {
            const int index = ib * nkb + ih;
            occ += (std::conj(becp[index]) * becp[index]).real();
        }
        Mi_z += sign * weight * occ;
    }

    // band0: |0.8|^2 + |0.6|^2 + |0.3+0.4i|^2 = 0.64 + 0.36 + 0.25 = 1.25, w=1.0
    // band1: |0.5|^2 + 0 + |1.0|^2 = 0.25 + 1.0 = 1.25, w=0.5
    // Mi_z = 1*1.25 + 0.5*1.25 = 1.875
    EXPECT_NEAR(Mi_z, 1.875, 1e-12);
}

TEST_F(DeltaSpinPwTest, MiPw_Npol2_PureZMag)
{
    // npol=2: construct becp so that only z-component is nonzero
    // becp_up = (a, 0), becp_dn = (0, 0)
    // occ[0] = |a|^2, occ[1]=0, occ[2]=0, occ[3]=0
    // Mi.z = w*(occ0-occ3) = w*|a|^2, Mi.x = 0, Mi.y = 0
    const int nkb = 1;
    const int nbands = 1;
    const double weight = 1.0;

    std::vector<std::complex<double>> becp(nbands * 2 * nkb, {0.0, 0.0});
    becp[0] = {0.7, 0.0};       // becp_up
    becp[0 + nkb] = {0.0, 0.0}; // becp_dn

    double Mi_x = 0.0, Mi_y = 0.0, Mi_z = 0.0;
    std::complex<double> occ[4] = {{0,0},{0,0},{0,0},{0,0}};
    occ[0] = std::conj(becp[0]) * becp[0];
    occ[1] = std::conj(becp[0]) * becp[0 + nkb];
    occ[2] = std::conj(becp[0 + nkb]) * becp[0];
    occ[3] = std::conj(becp[0 + nkb]) * becp[0 + nkb];

    Mi_z += weight * (occ[0] - occ[3]).real();
    Mi_x += weight * (occ[1] + occ[2]).real();
    Mi_y += weight * (occ[1] - occ[2]).imag();

    EXPECT_NEAR(Mi_z, 0.49, 1e-12);
    EXPECT_NEAR(Mi_x, 0.0, 1e-15);
    EXPECT_NEAR(Mi_y, 0.0, 1e-15);
}

TEST_F(DeltaSpinPwTest, MiPw_Npol2_PureXMag)
{
    // Construct becp so that only x-component is nonzero
    // becp_up = (a, 0), becp_dn = (a, 0) with same magnitude
    // occ[0] = |a|^2, occ[1] = |a|^2, occ[2] = |a|^2, occ[3] = |a|^2
    // Mi.z = w*(occ0-occ3) = 0
    // Mi.x = w*(occ1+occ2).real = w*2*|a|^2
    // Mi.y = w*(occ1-occ2).imag = 0
    const int nkb = 1;
    const int nbands = 1;
    const double weight = 1.0;
    const double a = 0.5;

    std::vector<std::complex<double>> becp(nbands * 2 * nkb, {0.0, 0.0});
    becp[0] = {a, 0.0};
    becp[0 + nkb] = {a, 0.0};

    std::complex<double> occ[4];
    occ[0] = std::conj(becp[0]) * becp[0];
    occ[1] = std::conj(becp[0]) * becp[0 + nkb];
    occ[2] = std::conj(becp[0 + nkb]) * becp[0];
    occ[3] = std::conj(becp[0 + nkb]) * becp[0 + nkb];

    double Mi_z = weight * (occ[0] - occ[3]).real();
    double Mi_x = weight * (occ[1] + occ[2]).real();
    double Mi_y = weight * (occ[1] - occ[2]).imag();

    EXPECT_NEAR(Mi_z, 0.0, 1e-15);
    EXPECT_NEAR(Mi_x, 0.5, 1e-12); // 2*0.25
    EXPECT_NEAR(Mi_y, 0.0, 1e-15);
}

TEST_F(DeltaSpinPwTest, MiPw_Npol2_MixedMag)
{
    // General becp: verify all three components
    const int nkb = 1;
    const int nbands = 1;
    const double weight = 1.0;

    std::vector<std::complex<double>> becp(nbands * 2 * nkb, {0.0, 0.0});
    becp[0] = {0.8, 0.0};        // becp_up
    becp[0 + nkb] = {0.0, 0.6};  // becp_dn

    std::complex<double> occ[4];
    occ[0] = std::conj(becp[0]) * becp[0];           // 0.64
    occ[1] = std::conj(becp[0]) * becp[0 + nkb];     // 0.8*(0,0.6) = (0, 0.48)
    occ[2] = std::conj(becp[0 + nkb]) * becp[0];     // (0,-0.6)*0.8 = (0, -0.48)
    occ[3] = std::conj(becp[0 + nkb]) * becp[0 + nkb]; // 0.36

    double Mi_z = weight * (occ[0] - occ[3]).real();
    double Mi_x = weight * (occ[1] + occ[2]).real();
    double Mi_y = weight * (occ[1] - occ[2]).imag();

    EXPECT_NEAR(Mi_z, 0.28, 1e-12);  // 0.64 - 0.36
    EXPECT_NEAR(Mi_x, 0.0, 1e-15);   // (0,0.48)+(0,-0.48) = 0
    EXPECT_NEAR(Mi_y, 0.96, 1e-12);  // imag((0,0.48)-(0,-0.48)) = imag(0,0.96) = 0.96
}

TEST_F(DeltaSpinPwTest, MiPw_MultiAtom_BeginIhOffset)
{
    // Two atoms with different nh, verify begin_ih offset
    const int nat = 2;
    const int nh_0 = 2, nh_1 = 1;
    const int nkb = nh_0 + nh_1; // 3
    const int nbands = 1;
    const double weight = 1.0;
    const int sign = 1;

    std::vector<std::complex<double>> becp(nbands * nkb, {0.0, 0.0});
    // atom 0: ih=0,1
    becp[0] = {1.0, 0.0}; // |becp|^2 = 1.0
    becp[1] = {0.0, 1.0}; // |becp|^2 = 1.0
    // atom 1: ih=2
    becp[2] = {0.5, 0.5}; // |becp|^2 = 0.5

    int nh_iat[2] = {nh_0, nh_1};
    double Mi_z[2] = {0.0, 0.0};

    for(int ib = 0; ib < nbands; ib++)
    {
        int begin_ih = 0;
        for(int iat = 0; iat < nat; iat++)
        {
            double occ = 0.0;
            for(int ih = 0; ih < nh_iat[iat]; ih++)
            {
                const int index = ib * nkb + begin_ih + ih;
                occ += (std::conj(becp[index]) * becp[index]).real();
            }
            Mi_z[iat] += sign * weight * occ;
            begin_ih += nh_iat[iat];
        }
    }

    EXPECT_NEAR(Mi_z[0], 2.0, 1e-12); // 1.0 + 1.0
    EXPECT_NEAR(Mi_z[1], 0.5, 1e-12); // 0.5
}

// =====================================================================
// cal_mw_from_lambda: magnetization re-accumulation from becp_tmp
// =====================================================================

TEST_F(DeltaSpinPwTest, MwFromLambda_Npol2_Accumulation)
{
    // Same formula as cal_Mi_pw npol=2, but from becp_tmp
    const int nkb = 1;
    const int nbands = 1;
    const int npol = 2;
    const int nk = 2;
    const double weights[2] = {1.0, 0.5};

    const int size_becp = nbands * nkb * npol;
    std::vector<std::complex<double>> becp_tmp(size_becp * nk, {0.0, 0.0});
    // k=0
    becp_tmp[0] = {0.8, 0.0};       // becp_up
    becp_tmp[0 + nkb] = {0.0, 0.6}; // becp_dn
    // k=1
    becp_tmp[size_becp + 0] = {0.6, 0.0};
    becp_tmp[size_becp + 0 + nkb] = {0.0, 0.8};

    double Mi_x = 0.0, Mi_y = 0.0, Mi_z = 0.0;
    int nh_iat[1] = {1};

    for(int ik = 0; ik < nk; ik++)
    {
        const std::complex<double>* becp = &becp_tmp[ik * size_becp];
        for(int ib = 0; ib < nbands; ib++)
        {
            const double weight = weights[ik];
            int begin_ih = 0;
            for(int iat = 0; iat < 1; iat++)
            {
                std::complex<double> occ[4] = {{0,0},{0,0},{0,0},{0,0}};
                for(int ih = 0; ih < nh_iat[iat]; ih++)
                {
                    const int index = ib * npol * nkb + begin_ih + ih;
                    occ[0] += std::conj(becp[index]) * becp[index];
                    occ[1] += std::conj(becp[index]) * becp[index + nkb];
                    occ[2] += std::conj(becp[index + nkb]) * becp[index];
                    occ[3] += std::conj(becp[index + nkb]) * becp[index + nkb];
                }
                Mi_x += weight * (occ[1] + occ[2]).real();
                Mi_y += weight * (occ[1] - occ[2]).imag();
                Mi_z += weight * (occ[0] - occ[3]).real();
                begin_ih += nh_iat[iat];
            }
        }
    }

    // k=0, w=1.0: occ0=0.64, occ3=0.36 => dz=0.28, occ1=(0,0.48), occ2=(0,-0.48) => dx=0, dy=0.96
    // k=1, w=0.5: occ0=0.36, occ3=0.64 => dz=-0.28*0.5=-0.14, occ1=(0,0.48), occ2=(0,-0.48) => dy=0.96*0.5=0.48
    EXPECT_NEAR(Mi_z, 0.14, 1e-12);  // 0.28 - 0.14
    EXPECT_NEAR(Mi_x, 0.0, 1e-15);
    EXPECT_NEAR(Mi_y, 1.44, 1e-12);  // 0.96 + 0.48
}

// =====================================================================
// DeltaHcc gemm contribution: h_tmp += becp^H * ps
// =====================================================================

TEST_F(DeltaSpinPwTest, DeltaHcc_GemmContribution)
{
    // Verify h_tmp += becp^H * ps for a small 2x2 case
    // becp: (npm x nbands), ps: (npm x nbands)
    // h_tmp += becp^H * ps = (nbands x npm) * (npm x nbands)
    const int nbands = 2;
    const int npm = 2; // nkb * npol

    // becp^H means conjugate transpose
    std::vector<std::complex<double>> becp = {
        {1.0, 0.0}, {0.0, 1.0},  // column 0: becp[0,0], becp[1,0]
        {0.5, 0.0}, {0.0, -0.5}  // column 1: becp[0,1], becp[1,1]
    };
    std::vector<std::complex<double>> ps = {
        {0.5, 0.0}, {0.0, 0.5},
        {0.3, 0.0}, {0.0, -0.3}
    };

    // Manual: h_tmp[i,j] += sum_k conj(becp[k,i]) * ps[k,j]
    // becp stored column-major as becp[k + i*npm], ps stored as ps[k + j*npm]
    std::vector<std::complex<double>> h_tmp(nbands * nbands, {0.0, 0.0});
    for(int i = 0; i < nbands; i++)
    {
        for(int j = 0; j < nbands; j++)
        {
            for(int k = 0; k < npm; k++)
            {
                h_tmp[i * nbands + j] += std::conj(becp[k + i * npm]) * ps[k + j * npm];
            }
        }
    }

    // h[0,0] = conj(1)*0.5 + conj(0,1)*(0,0.5) = 0.5 + (0,-1)*(0,0.5) = 0.5 + 0.5 = 1.0
    EXPECT_NEAR(h_tmp[0].real(), 1.0, 1e-12);
    EXPECT_NEAR(h_tmp[0].imag(), 0.0, 1e-12);

    // h[0,1] = conj(1)*0.3 + conj(0,1)*(0,-0.3) = 0.3 + (0,-1)*(0,-0.3) = 0.3 + (-0.3) = 0
    EXPECT_NEAR(h_tmp[1].real(), 0.0, 1e-12);
    EXPECT_NEAR(h_tmp[1].imag(), 0.0, 1e-12);

    // h[1,0] = conj(0.5)*0.5 + conj(0,-0.5)*(0,0.5) = 0.25 + (0,0.5)*(0,0.5) = 0.25 + (-0.25) = 0
    EXPECT_NEAR(h_tmp[2].real(), 0.0, 1e-12);
    EXPECT_NEAR(h_tmp[2].imag(), 0.0, 1e-12);

    // h[1,1] = conj(0.5)*0.3 + conj(0,-0.5)*(0,-0.3) = 0.15 + (0,0.5)*(0,-0.3) = 0.15 + 0.15 = 0.3
    EXPECT_NEAR(h_tmp[3].real(), 0.3, 1e-12);
    EXPECT_NEAR(h_tmp[3].imag(), 0.0, 1e-12);
}
