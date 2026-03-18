#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "../module_ewald/H_Ewald_pw.h"
#include "../module_ewald/dnrm2.h"
#include "source_base/matrix3.h"
#include <vector>

/************************************************
 *  unit test of H_Ewald_pw::rgen
 ***********************************************/

/**
 * - Tested Functions:
 *   - H_Ewald_pw::rgen():
 *      - Generates lattice vectors R such that |R - dtau| <= rmax,
 *        and returns them sorted by ascending magnitude.
 *      - Tested cases:
 *        1. rmax == 0.0: no vectors returned.
 *        2. Simple cubic cell, small rmax: correct count + sorted order.
 *        3. Large rmax exceeding original fixed mxr=200 limit: verifies
 *           that the dynamic mxr sizing introduced in the bug fix works
 *           correctly and does not overflow the allocated arrays.
 */

class RgenTest : public ::testing::Test
{
  protected:
    // Simple cubic unit cell: latvec = G = identity
    ModuleBase::Matrix3 latvec;
    ModuleBase::Matrix3 G;
    ModuleBase::Vector3<double> dtau;

    void SetUp() override
    {
        // Unit cubic cell: direct and reciprocal lattice vectors are identity
        latvec = ModuleBase::Matrix3(1, 0, 0, 0, 1, 0, 0, 0, 1);
        G = ModuleBase::Matrix3(1, 0, 0, 0, 1, 0, 0, 0, 1);
        dtau = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
    }
};

TEST_F(RgenTest, ZeroRmax)
{
    // When rmax==0 the function should return immediately with nrm=0
    const int mxr_test = 10;
    std::vector<ModuleBase::Vector3<double>> r(mxr_test);
    std::vector<double> r2(mxr_test);
    std::vector<int> irr(mxr_test);
    int nrm = 0;

    H_Ewald_pw::rgen(dtau, 0.0, irr.data(), latvec, G, r.data(), r2.data(), mxr_test, nrm);

    EXPECT_EQ(nrm, 0);
}

TEST_F(RgenTest, SimpleCubicNearestNeighbors)
{
    // rmax = 1.5 captures nearest (d=1) and next-nearest (d=sqrt(2)~1.414)
    // neighbors: 6 + 12 = 18 vectors total.
    const double rmax = 1.5;
    const int mxr_test = 50;
    std::vector<ModuleBase::Vector3<double>> r(mxr_test);
    std::vector<double> r2(mxr_test);
    std::vector<int> irr(mxr_test);
    int nrm = 0;

    H_Ewald_pw::rgen(dtau, rmax, irr.data(), latvec, G, r.data(), r2.data(), mxr_test, nrm);

    EXPECT_EQ(nrm, 18);

    // Vectors must be sorted in ascending order of |r|^2
    for (int i = 1; i < nrm; ++i)
    {
        EXPECT_LE(r2[i - 1], r2[i]);
    }

    // All returned vectors must lie strictly inside the sphere
    for (int i = 0; i < nrm; ++i)
    {
        EXPECT_LE(r2[i], rmax * rmax + 1.0e-10);
        EXPECT_GT(r2[i], 1.0e-10);
    }
}

TEST_F(RgenTest, SimpleCubicNonZeroDtau)
{
    // rgen computes t = R - dtau for each lattice vector R=(i,j,k)*latvec,
    // and excludes vectors with |t|^2 < 1e-10 (i.e. R == dtau).
    // With dtau=(0.5,0,0) and rmax=0.6, two vectors qualify:
    //   R=(0,0,0): t = (0,0,0)-(0.5,0,0) = (-0.5,0,0), |t|^2=0.25 <= 0.36
    //   R=(1,0,0): t = (1,0,0)-(0.5,0,0) = ( 0.5,0,0), |t|^2=0.25 <= 0.36
    // No lattice point coincides with dtau, so neither is excluded.
    const double rmax = 0.6;
    const int mxr_test = 10;
    dtau = ModuleBase::Vector3<double>(0.5, 0.0, 0.0);
    std::vector<ModuleBase::Vector3<double>> r(mxr_test);
    std::vector<double> r2(mxr_test);
    std::vector<int> irr(mxr_test);
    int nrm = 0;

    H_Ewald_pw::rgen(dtau, rmax, irr.data(), latvec, G, r.data(), r2.data(), mxr_test, nrm);

    EXPECT_EQ(nrm, 2);
    for (int i = 0; i < nrm; ++i)
    {
        EXPECT_NEAR(r2[i], 0.25, 1.0e-10);
    }
}

TEST_F(RgenTest, LargeRmaxExceedsOriginalLimit)
{
    // rmax=4.0 on a unit cubic cell yields ~499 r-vectors, well above the
    // old fixed limit of mxr=200 that caused the buffer overflow.
    // This test verifies that with a properly sized mxr the function
    // completes without error.
    const double rmax = 4.0;

    // Replicate the dynamic mxr computation from compute_ewald()
    double bg1[3];
    bg1[0] = G.e11; bg1[1] = G.e12; bg1[2] = G.e13;
    int nm1 = (int)(dnrm2(3, bg1, 1) * rmax + 2);
    bg1[0] = G.e21; bg1[1] = G.e22; bg1[2] = G.e23;
    int nm2 = (int)(dnrm2(3, bg1, 1) * rmax + 2);
    bg1[0] = G.e31; bg1[1] = G.e32; bg1[2] = G.e33;
    int nm3 = (int)(dnrm2(3, bg1, 1) * rmax + 2);
    const int mxr_test = (2 * nm1 + 1) * (2 * nm2 + 1) * (2 * nm3 + 1);

    std::vector<ModuleBase::Vector3<double>> r(mxr_test);
    std::vector<double> r2(mxr_test);
    std::vector<int> irr(mxr_test);
    int nrm = 0;

    H_Ewald_pw::rgen(dtau, rmax, irr.data(), latvec, G, r.data(), r2.data(), mxr_test, nrm);

    // Must exceed the old hard-coded limit that caused the crash
    EXPECT_GT(nrm, 200);

    // All returned vectors lie within the sphere
    for (int i = 0; i < nrm; ++i)
    {
        EXPECT_LE(r2[i], rmax * rmax + 1.0e-10);
        EXPECT_GT(r2[i], 1.0e-10);
    }

    // Vectors are sorted in ascending order of |r|^2
    for (int i = 1; i < nrm; ++i)
    {
        EXPECT_LE(r2[i - 1], r2[i]);
    }
}
