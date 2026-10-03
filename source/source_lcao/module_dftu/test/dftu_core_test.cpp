#include "gtest/gtest.h"
#include <cmath>
#include <complex>
#include <vector>
#include <numeric>
#include <algorithm>

/***********************************************************************
 * Unit tests for DFT+U core algorithms.
 *
 * These tests target the most complex and bug-prone logic:
 * 1. pot_onsite effective potential formula (UForm::dud_fll, FLL)
 * 2. Energy correction and double-counting terms
 ***********************************************************************/

// =====================================================================
// 1. pot_onsite effective potential formula (UForm::dud_fll, FLL)
//
// pot_onsite[m0,m1] = U * (0.5*delta(m0,m1) - occ_mat[m0,m1])  (diagonal)
// pot_onsite[m0,m1] = -U * occ_mat[m0,m1]                       (off-diagonal)
// =====================================================================

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

// =====================================================================
// 4. Energy correction formula
//
// E_U = 0.5 * U * sum_spin [Tr(n) - Tr(n^2)]
// =====================================================================

class EnergyCorrectionTest : public ::testing::Test
{
  protected:
    static double compute_energy(const std::vector<double>& occ_mat_flat, int m_size, double U)
    {
        double nm_trace = 0.0, nm2_trace = 0.0;
        for (int m0 = 0; m0 < m_size; m0++)
        {
            nm_trace += occ_mat_flat[m0 * m_size + m0];
            for (int m1 = 0; m1 < m_size; m1++)
                nm2_trace += occ_mat_flat[m0 * m_size + m1] * occ_mat_flat[m1 * m_size + m0];
        }
        return 0.5 * U * (nm_trace - nm2_trace);
    }
};

TEST_F(EnergyCorrectionTest, HalfFilled_DOrbital)
{
    const int m_size = 5;
    std::vector<double> occ_mat(m_size * m_size, 0.0);
    for (int m = 0; m < m_size; m++)
        occ_mat[m * m_size + m] = 0.5;

    double energy = compute_energy(occ_mat, m_size, 4.0);
    // Tr(n) = 2.5, Tr(n^2) = 1.25, E = 0.5 * 4 * 1.25 = 2.5
    EXPECT_DOUBLE_EQ(energy, 2.5);
}

TEST_F(EnergyCorrectionTest, OffDiagonal_Contribution)
{
    const int m_size = 2;
    std::vector<double> occ_mat = {
        0.3, 0.1,
        0.1, 0.3
    };

    double energy = compute_energy(occ_mat, m_size, 4.0);
    // Tr(n) = 0.6, Tr(n^2) = 0.3^2 + 0.1^2 + 0.1^2 + 0.3^2 = 0.20
    // E = 0.5 * 4 * (0.6 - 0.20) = 0.8
    EXPECT_DOUBLE_EQ(energy, 0.8);
}

TEST_F(EnergyCorrectionTest, DoubleCounting_Energy)
{
    // E_dc = sum_{m1,m2,spin} pot_onsite[m1,m2] * n[m2,m1]
    const int m_size = 3;
    double U = 4.0;
    std::vector<double> occ_mat = {
        0.5, 0.0, 0.0,
        0.0, 0.3, 0.0,
        0.0, 0.0, 0.2
    };

    double e_dc = 0.0;
    for (int m1 = 0; m1 < m_size; m1++)
        for (int m2 = 0; m2 < m_size; m2++)
        {
            double pot_onsite = (m1 == m2) ? U * (0.5 - occ_mat[m1 * m_size + m2])
                                   : -U * occ_mat[m1 * m_size + m2];
            e_dc += pot_onsite * occ_mat[m2 * m_size + m1];
        }

    // Only diagonal: m=0: 0*0.5=0, m=1: 0.8*0.3=0.24, m=2: 1.2*0.2=0.24
    EXPECT_NEAR(e_dc, 0.48, 1e-14);
}
