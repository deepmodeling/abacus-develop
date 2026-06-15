#include "gtest/gtest.h"
#include <cmath>
#include <complex>
#include <vector>
#include <cstring>

#include "source_pw/module_pwdft/kernels/onsite_op.h"
#include "source_base/module_device/types.h"

/***********************************************************************
 * Unit tests for onsite_ps_op kernel (CPU implementation).
 *
 * Tests the actual kernel functor for both DFT+U and DeltaSpin paths,
 * covering npol=1 (collinear) and npol=2 (non-collinear) branches.
 ***********************************************************************/

using complexd = std::complex<double>;

// =====================================================================
// 1. DeltaSpin kernel (npol=1 branch)
//
// For npol=1 (nspin=1 or nspin=2 collinear):
//   ps[ip * npm + ib] += lambda_array[iat] * becp[ib * tnp + ip]
// =====================================================================

class OnsitePsDeltaSpinNpol1Test : public ::testing::Test
{
  protected:
    hamilt::onsite_ps_op<double, base_device::DEVICE_CPU> kernel;

    void SetUp() override {}
};

TEST_F(OnsitePsDeltaSpinNpol1Test, SingleBandSingleProj)
{
    const int npm = 1, npol = 1, tnp = 1;

    std::vector<int> ip_iat = {0};
    std::vector<complexd> lambda_array = {2.0}; // lambda for atom 0
    std::vector<complexd> becp = {complexd(1.0, 0.5)};
    std::vector<complexd> ps = {0.0};

    kernel(nullptr, npm, npol, ip_iat.data(), tnp,
           lambda_array.data(), ps.data(), becp.data());

    // ps[0] += 2.0 * (complexd(1.0, 0.5)) = complexd(2.0, 1.0)
    EXPECT_NEAR(ps[0].real(), 2.0, 1e-15);
    EXPECT_NEAR(ps[0].imag(), 1.0, 1e-15);
}

TEST_F(OnsitePsDeltaSpinNpol1Test, MultiBandMultiProj)
{
    const int npm = 2, npol = 1, tnp = 2;

    std::vector<int> ip_iat = {0, 1};
    std::vector<complexd> lambda_array = {3.0, 4.0}; // lambda for atom 0 and 1
    // becp[ib * tnp + ip]: band 0 has proj 0,1; band 1 has proj 0,1
    std::vector<complexd> becp = {
        complexd(1.0, 0.0), complexd(0.5, 0.0),  // band 0: proj 0, 1
        0.0 + complexd(0.0, 1.0), 0.0 + complexd(0.0, 2.0)   // band 1: proj 0, 1
    };
    std::vector<complexd> ps = {0.0, 0.0, 0.0, 0.0}; // tnp * npm = 4

    kernel(nullptr, npm, npol, ip_iat.data(), tnp,
           lambda_array.data(), ps.data(), becp.data());

    // ps[ip * npm + ib] += lambda[iat] * becp[ib * tnp + ip]
    // ib=0,ip=0: ps[0] += 3.0 * 1.0 = 3.0
    // ib=0,ip=1: ps[2] += 4.0 * 0.5 = 2.0
    // ib=1,ip=0: ps[1] += 3.0 * (0+1i) = 0+3i
    // ib=1,ip=1: ps[3] += 4.0 * (0+2i) = 0+8i
    EXPECT_NEAR(ps[0].real(), 3.0, 1e-15);
    EXPECT_NEAR(ps[1].imag(), 3.0, 1e-15);
    EXPECT_NEAR(ps[2].real(), 2.0, 1e-15);
    EXPECT_NEAR(ps[3].imag(), 8.0, 1e-15);
}

TEST_F(OnsitePsDeltaSpinNpol1Test, AccumulatesOnExistingPs)
{
    // Kernel uses +=, so it should accumulate on existing values
    const int npm = 1, npol = 1, tnp = 1;

    std::vector<int> ip_iat = {0};
    std::vector<complexd> lambda_array = {2.0};
    std::vector<complexd> becp = {1.0};
    std::vector<complexd> ps = {5.0}; // pre-existing value

    kernel(nullptr, npm, npol, ip_iat.data(), tnp,
           lambda_array.data(), ps.data(), becp.data());

    // ps[0] = 5.0 + 2.0 * 1.0 = 7.0
    EXPECT_NEAR(ps[0].real(), 7.0, 1e-15);
}

// =====================================================================
// 2. DeltaSpin kernel (npol=2 branch)
//
// For npol=2 (non-collinear):
//   ps[ip * npm + ib2] += lambda[4*iat+0] * becp[ib2*tnp + ip]
//                       + lambda[4*iat+2] * becp[ib2*tnp + ip + tnp]
//   ps[ip * npm + ib2+1] += lambda[4*iat+1] * becp[ib2*tnp + ip]
//                         + lambda[4*iat+3] * becp[ib2*tnp + ip + tnp]
// =====================================================================

class OnsitePsDeltaSpinNpol2Test : public ::testing::Test
{
  protected:
    hamilt::onsite_ps_op<double, base_device::DEVICE_CPU> kernel;

    void SetUp() override {}
};

TEST_F(OnsitePsDeltaSpinNpol2Test, SingleBandSingleProj)
{
    const int npm = 2, npol = 2, tnp = 1; // npm/npol = 1 band

    std::vector<int> ip_iat = {0};
    // lambda coefficients for atom 0: Pauli matrix elements
    // [lambda_z, lambda_x+iy, lambda_x-iy, -lambda_z]
    std::vector<complexd> lambda_array = {
        3.0, complexd(1.0, 2.0), complexd(1.0, -2.0), -3.0
    };
    // becp: 2 rows (spin up/down) x 1 proj
    std::vector<complexd> becp = {
        1.0 + complexd(0.0, 0.0),  // spin up
        complexd(0.0, 1.0)   // spin down
    };
    std::vector<complexd> ps = {0.0, 0.0};

    kernel(nullptr, npm, npol, ip_iat.data(), tnp,
           lambda_array.data(), ps.data(), becp.data());

    // ps[0] += lambda[0] * becp[0] + lambda[2] * becp[1]
    //       = 3.0 * 1.0 + (1-complexd(0.0, 2.0)) * complexd(0.0, 1.0)
    //       = 3.0 + (complexd(0.0, 1.0) + 2) = 5.0 + complexd(0.0, 1.0)
    // ps[1] += lambda[1] * becp[0] + lambda[3] * becp[1]
    //       = (1+complexd(0.0, 2.0)) * 1.0 + (-3.0) * complexd(0.0, 1.0)
    //       = 1+complexd(0.0, 2.0) - complexd(0.0, 3.0) = 1.0 - complexd(0.0, 1.0)
    EXPECT_NEAR(ps[0].real(), 5.0, 1e-15);
    EXPECT_NEAR(ps[0].imag(), 1.0, 1e-15);
    EXPECT_NEAR(ps[1].real(), 1.0, 1e-15);
    EXPECT_NEAR(ps[1].imag(), -1.0, 1e-15);
}

TEST_F(OnsitePsDeltaSpinNpol2Test, MultiBand)
{
    const int npm = 4, npol = 2, tnp = 1; // npm/npol = 2 bands

    std::vector<int> ip_iat = {0};
    std::vector<complexd> lambda_array = {2.0, 0.0, 0.0, -2.0};
    std::vector<complexd> becp = {
        1.0, 0.0,    // band 0: up, down
        0.0, 1.0     // band 1: up, down
    };
    std::vector<complexd> ps = {0.0, 0.0, 0.0, 0.0};

    kernel(nullptr, npm, npol, ip_iat.data(), tnp,
           lambda_array.data(), ps.data(), becp.data());

    // Band 0 (ib=0, ib2=0):
    //   ps[0] += 2.0*1.0 + 0.0*0.0 = 2.0
    //   ps[1] += 0.0*1.0 + (-2.0)*0.0 = 0.0
    // Band 1 (ib=1, ib2=2):
    //   ps[2] += 2.0*0.0 + 0.0*1.0 = 0.0
    //   ps[3] += 0.0*0.0 + (-2.0)*1.0 = -2.0
    EXPECT_NEAR(ps[0].real(), 2.0, 1e-15);
    EXPECT_NEAR(ps[3].real(), -2.0, 1e-15);
}

// =====================================================================
// 3. DFT+U kernel (npol=1 branch)
//
// For npol=1:
//   for each ip:
//     m1 = ip_m[ip], if m1 < 0 continue
//     iat = ip_iat[ip], vu_iat = vu + vu_begin_iat[iat]
//     tlp1 = 2*orb_l + 1
//     for ip2 in [ip-m1, ip-m1+tlp1):
//       m2 = ip_m[ip2]
//       ps[ip * npm + ib] += vu_iat[m1*tlp1 + m2] * becp[ib * tnp + ip2]
// =====================================================================

class OnsitePsDftuNpol1Test : public ::testing::Test
{
  protected:
    hamilt::onsite_ps_op<double, base_device::DEVICE_CPU> kernel;

    void SetUp() override {}
};

TEST_F(OnsitePsDftuNpol1Test, SingleBandSingleAtom_DOrbital)
{
    // 1 band, 1 atom with d-orbital (5 projectors), all correlated
    const int npm = 1, npol = 1, tnp = 5;

    std::vector<int> orb_l_iat = {2}; // d-orbital (l=2)
    std::vector<int> ip_iat = {0, 0, 0, 0, 0}; // all belong to atom 0
    std::vector<int> ip_m = {0, 1, 2, 3, 4}; // m indices
    std::vector<int> vu_begin_iat = {0}; // VU starts at index 0

    // VU matrix for d-orbital (5x5), row-major
    std::vector<complexd> vu(25, 0.0);
    vu[0] = 1.0; // VU[0,0]
    vu[6] = 2.0; // VU[1,1]
    vu[12] = 3.0; // VU[2,2]

    std::vector<complexd> becp = {1.0, 0.5, 0.3, 0.2, 0.1}; // 5 projectors
    std::vector<complexd> ps(5, 0.0); // tnp * npm = 5

    kernel(nullptr, npm, npol,
           orb_l_iat.data(), ip_iat.data(), ip_m.data(), vu_begin_iat.data(),
           tnp, vu.data(), ps.data(), becp.data());

    // For ip=0 (m1=0): ip2 ranges from 0 to 5
    //   ps[0] += VU[0,0]*becp[0] + VU[0,1]*becp[1] + ...
    //          = 1.0*1.0 + 0 + 0 + 0 + 0 = 1.0
    // For ip=1 (m1=1): ip2 ranges from 0 to 5
    //   ps[1] += VU[1,0]*becp[0] + VU[1,1]*becp[1] + ...
    //          = 0 + 2.0*0.5 + 0 + 0 + 0 = 1.0
    // For ip=2 (m1=2): ps[2] += 3.0*0.3 = 0.9
    EXPECT_NEAR(ps[0].real(), 1.0, 1e-15);
    EXPECT_NEAR(ps[1].real(), 1.0, 1e-15);
    EXPECT_NEAR(ps[2].real(), 0.9, 1e-15);
    EXPECT_NEAR(ps[3].real(), 0.0, 1e-15);
    EXPECT_NEAR(ps[4].real(), 0.0, 1e-15);
}

TEST_F(OnsitePsDftuNpol1Test, OffDiagonalVU)
{
    // Test off-diagonal VU elements
    const int npm = 1, npol = 1, tnp = 3; // p-orbital

    std::vector<int> orb_l_iat = {1}; // p-orbital
    std::vector<int> ip_iat = {0, 0, 0};
    std::vector<int> ip_m = {0, 1, 2};
    std::vector<int> vu_begin_iat = {0};

    // VU with off-diagonal: VU[0,1] = 0.5, VU[1,0] = 0.5
    std::vector<complexd> vu(9, 0.0);
    vu[1] = 0.5; // VU[0,1]
    vu[3] = 0.5; // VU[1,0]

    std::vector<complexd> becp = {1.0, 2.0, 3.0};
    std::vector<complexd> ps(3, 0.0);

    kernel(nullptr, npm, npol,
           orb_l_iat.data(), ip_iat.data(), ip_m.data(), vu_begin_iat.data(),
           tnp, vu.data(), ps.data(), becp.data());

    // ip=0 (m1=0): ip2 from 0 to 3
    //   ps[0] += VU[0,0]*becp[0] + VU[0,1]*becp[1] + VU[0,2]*becp[2]
    //          = 0*1.0 + 0.5*2.0 + 0*3.0 = 1.0
    // ip=1 (m1=1):
    //   ps[1] += VU[1,0]*becp[0] + VU[1,1]*becp[1] + VU[1,2]*becp[2]
    //          = 0.5*1.0 + 0*2.0 + 0*3.0 = 0.5
    EXPECT_NEAR(ps[0].real(), 1.0, 1e-15);
    EXPECT_NEAR(ps[1].real(), 0.5, 1e-15);
}

TEST_F(OnsitePsDftuNpol1Test, NonCorrelatedProjector_MMinus1)
{
    // Test that ip_m = -1 projectors are skipped
    const int npm = 1, npol = 1, tnp = 4;

    std::vector<int> orb_l_iat = {1}; // p-orbital
    std::vector<int> ip_iat = {0, 0, 0, 0};
    // First projector is not correlated (m=-1), rest are p-type (m=0,1,2)
    std::vector<int> ip_m = {-1, 0, 1, 2};
    std::vector<int> vu_begin_iat = {0};

    std::vector<complexd> vu(9, 1.0); // all VU = 1.0
    std::vector<complexd> becp = {1.0, 1.0, 1.0, 1.0};
    std::vector<complexd> ps(4, 0.0);

    kernel(nullptr, npm, npol,
           orb_l_iat.data(), ip_iat.data(), ip_m.data(), vu_begin_iat.data(),
           tnp, vu.data(), ps.data(), becp.data());

    // ip=0: m1=-1, skipped
    // ip=1,2,3: should have contributions
    EXPECT_NEAR(ps[0].real(), 0.0, 1e-15); // skipped
    EXPECT_NEAR(ps[1].real(), 3.0, 1e-15); // sum of VU[1,*]*becp
    EXPECT_NEAR(ps[2].real(), 3.0, 1e-15);
    EXPECT_NEAR(ps[3].real(), 3.0, 1e-15);
}

TEST_F(OnsitePsDftuNpol1Test, MultiBand_DOrbital)
{
    const int npm = 2, npol = 1, tnp = 5;

    std::vector<int> orb_l_iat = {2};
    std::vector<int> ip_iat = {0, 0, 0, 0, 0};
    std::vector<int> ip_m = {0, 1, 2, 3, 4};
    std::vector<int> vu_begin_iat = {0};

    std::vector<complexd> vu(25, 0.0);
    vu[0] = 2.0; // VU[0,0]

    // Band 0: becp[0..4], Band 1: becp[5..9]
    std::vector<complexd> becp = {
        1.0, 0.0, 0.0, 0.0, 0.0,  // band 0
        0.0, 1.0, 0.0, 0.0, 0.0   // band 1
    };
    std::vector<complexd> ps(10, 0.0); // tnp * npm = 10

    kernel(nullptr, npm, npol,
           orb_l_iat.data(), ip_iat.data(), ip_m.data(), vu_begin_iat.data(),
           tnp, vu.data(), ps.data(), becp.data());

    // Band 0: ps[0] += VU[0,0] * becp[0] = 2.0 * 1.0 = 2.0
    // Band 1: ps[5] += VU[0,0] * becp[5] = 2.0 * 0.0 = 0.0
    // Wait, becp indexing: becp[ib * tnp + ip2]
    // For band 1, ib=1: becp[1*5 + 0] = becp[5] = 0.0
    EXPECT_NEAR(ps[0].real(), 2.0, 1e-15);
    EXPECT_NEAR(ps[5].real(), 0.0, 1e-15);
}

// =====================================================================
// 4. DFT+U kernel (npol=2 branch)
//
// For npol=2:
//   ps[ip * npm + ib2] += vu_iat[index_mm] * becp[ib2*tnp + ip2]
//                       + vu_iat[index_mm + 2*tlp1^2] * becp[ib2*tnp + ip2 + tnp]
//   ps[ip * npm + ib2+1] += vu_iat[index_mm + tlp1^2] * becp[ib2*tnp + ip2]
//                         + vu_iat[index_mm + 3*tlp1^2] * becp[ib2*tnp + ip2 + tnp]
//
// where index_mm = m1 * tlp1 + m2
// =====================================================================

class OnsitePsDftuNpol2Test : public ::testing::Test
{
  protected:
    hamilt::onsite_ps_op<double, base_device::DEVICE_CPU> kernel;

    void SetUp() override {}
};

TEST_F(OnsitePsDftuNpol2Test, SingleBandSingleAtom_Porbital)
{
    // 1 band pair (npm/npol = 1), p-orbital (3 projectors)
    const int npm = 2, npol = 2, tnp = 3;

    std::vector<int> orb_l_iat = {1}; // p-orbital
    std::vector<int> ip_iat = {0, 0, 0};
    std::vector<int> ip_m = {0, 1, 2};
    std::vector<int> vu_begin_iat = {0};

    // VU: 4 blocks of 3x3 = 36 elements
    std::vector<complexd> vu(36, 0.0);
    // Block 0 (Pauli I): VU[0,0] = 2.0
    vu[0] = 2.0;
    // Block 1 (Pauli X): VU[0,0] = 1.0
    vu[9] = 1.0; // tlp1^2 = 9

    // becp: 2 rows (spin up/down) x 3 projectors
    std::vector<complexd> becp = {
        1.0, 0.0, 0.0,  // spin up
        0.0, 1.0, 0.0   // spin down
    };
    std::vector<complexd> ps = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // tnp * npm = 6

    kernel(nullptr, npm, npol,
           orb_l_iat.data(), ip_iat.data(), ip_m.data(), vu_begin_iat.data(),
           tnp, vu.data(), ps.data(), becp.data());

    // For ip=0 (m1=0): ip2 from 0 to 3
    //   ps[0] += vu[0]*becp[0] + vu[18]*becp[3] = 2.0*1.0 + 0*0.0 = 2.0
    //   ps[1] += vu[9]*becp[0] + vu[27]*becp[3] = 1.0*1.0 + 0*0.0 = 1.0
    EXPECT_NEAR(ps[0].real(), 2.0, 1e-15);
    EXPECT_NEAR(ps[1].real(), 1.0, 1e-15);
}

TEST_F(OnsitePsDftuNpol2Test, MultiBand)
{
    const int npm = 4, npol = 2, tnp = 3; // 2 band pairs

    std::vector<int> orb_l_iat = {1};
    std::vector<int> ip_iat = {0, 0, 0};
    std::vector<int> ip_m = {0, 1, 2};
    std::vector<int> vu_begin_iat = {0};

    std::vector<complexd> vu(36, 0.0);
    vu[0] = 2.0; // Block 0, VU[0,0]

    // becp: 2 spin x 3 proj x 2 bands = 12 elements
    // Layout: [band0_up, band0_dn, band1_up, band1_dn]
    std::vector<complexd> becp = {
        1.0, 0.5, 0.3,  // band 0: up
        0.0, 1.0, 0.2,  // band 0: down
        0.2, 0.0, 0.0,  // band 1: up
        0.0, 0.3, 0.0   // band 1: down
    };
    std::vector<complexd> ps(12, 0.0);

    kernel(nullptr, npm, npol,
           orb_l_iat.data(), ip_iat.data(), ip_m.data(), vu_begin_iat.data(),
           tnp, vu.data(), ps.data(), becp.data());

    // Band pair 0 (ib=0, ib2=0):
    //   ps[0] += vu[0]*becp[0] + vu[18]*becp[3] = 2.0*1.0 + 0 = 2.0
    // Band pair 1 (ib=1, ib2=2):
    //   ps[2] += vu[0]*becp[6] + vu[18]*becp[9] = 2.0*0.2 + 0 = 0.4
    EXPECT_NEAR(ps[0].real(), 2.0, 1e-15);
    EXPECT_NEAR(ps[2].real(), 0.4, 1e-15);
}

// =====================================================================
// 5. Edge cases and boundary conditions
// =====================================================================

class OnsitePsEdgeCasesTest : public ::testing::Test
{
  protected:
    hamilt::onsite_ps_op<double, base_device::DEVICE_CPU> kernel;

    void SetUp() override {}
};

TEST_F(OnsitePsEdgeCasesTest, EmptyBecp_NpmZero)
{
    const int npm = 0, npol = 1, tnp = 1;

    std::vector<int> ip_iat = {0};
    std::vector<complexd> lambda_array = {1.0};
    std::vector<complexd> becp;
    std::vector<complexd> ps;

    // Should not crash with npm=0
    EXPECT_NO_THROW(kernel(nullptr, npm, npol, ip_iat.data(), tnp,
                           lambda_array.data(), ps.data(), becp.data()));
}

TEST_F(OnsitePsEdgeCasesTest, ZeroLambda)
{
    const int npm = 2, npol = 1, tnp = 2;

    std::vector<int> ip_iat = {0, 1};
    std::vector<complexd> lambda_array = {0.0, 0.0};
    std::vector<complexd> becp = {1.0, 2.0, 3.0, 4.0};
    std::vector<complexd> ps = {10.0, 20.0, 30.0, 40.0};

    kernel(nullptr, npm, npol, ip_iat.data(), tnp,
           lambda_array.data(), ps.data(), becp.data());

    // Zero lambda should not change ps (0 += 0)
    EXPECT_NEAR(ps[0].real(), 10.0, 1e-15);
    EXPECT_NEAR(ps[1].real(), 20.0, 1e-15);
}

TEST_F(OnsitePsEdgeCasesTest, ComplexVU)
{
    // DFT+U with complex VU elements
    const int npm = 1, npol = 1, tnp = 2;

    std::vector<int> orb_l_iat = {0}; // s-orbital, but 2 projectors
    std::vector<int> ip_iat = {0, 0};
    std::vector<int> ip_m = {0, -1}; // first correlated, second not
    std::vector<int> vu_begin_iat = {0};

    std::vector<complexd> vu = {complexd(1.0, 2.0)}; // 1x1 VU matrix
    std::vector<complexd> becp = {complexd(0.5, 0.5), 1.0};
    std::vector<complexd> ps = {0.0, 0.0};

    kernel(nullptr, npm, npol,
           orb_l_iat.data(), ip_iat.data(), ip_m.data(), vu_begin_iat.data(),
           tnp, vu.data(), ps.data(), becp.data());

    // ps[0] += (1+complexd(0.0, 2.0)) * (0.5+0.5i) = 0.5+0.5i + complexd(0.0, 1.0)-1 = -0.5+1.5i
    EXPECT_NEAR(ps[0].real(), -0.5, 1e-15);
    EXPECT_NEAR(ps[0].imag(), 1.5, 1e-15);
    EXPECT_NEAR(ps[1].real(), 0.0, 1e-15); // m=-1, skipped
}
