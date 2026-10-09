#include "gtest/gtest.h"
#include <cmath>
#include <vector>

/***********************************************************************
 * Unit tests for the DFT+U on-site potential formula (dftu_nao_pots.cpp).
 *
 * pot_onsite[m0,m1] = U * (0.5*delta(m0,m1) - occ_mat[m0,m1])  (diagonal)
 * pot_onsite[m0,m1] = -U * occ_mat[m0,m1]                       (off-diagonal)
 ***********************************************************************/

static double compute_pot_onsite(double U_val, int m0, int m1, double occ_mat_val)
{
    if (m0 == m1)
        return U_val * (0.5 - occ_mat_val);
    else
        return -U_val * occ_mat_val;
}

class PotOnsitePotentialTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
};

TEST_F(PotOnsitePotentialTest, Diagonal_HalfFilled)
{
    double U = 4.0;
    double occ_mat = 0.5; // half-filled
    double pot_onsite = compute_pot_onsite(U, 0, 0, occ_mat);
    EXPECT_DOUBLE_EQ(pot_onsite, 0.0); // U * (0.5 - 0.5) = 0
}

TEST_F(PotOnsitePotentialTest, Diagonal_FullyOccupied)
{
    double U = 4.0;
    double occ_mat = 1.0; // fully occupied
    double pot_onsite = compute_pot_onsite(U, 0, 0, occ_mat);
    EXPECT_DOUBLE_EQ(pot_onsite, -2.0); // U * (0.5 - 1.0) = -2.0
}

TEST_F(PotOnsitePotentialTest, OffDiagonal)
{
    double U = 5.0;
    double occ_mat = 0.3;
    double pot_onsite = compute_pot_onsite(U, 0, 1, occ_mat);
    EXPECT_DOUBLE_EQ(pot_onsite, -1.5); // -U * occ_mat = -1.5
}
