#include "gtest/gtest.h"
#include <cmath>
#include <vector>

/***********************************************************************
 * Unit tests for the DFT+U energy correction formula (dftu_nao_energy.cpp).
 *
 * E_U = 0.5 * U * sum_spin [Tr(n) - Tr(n^2)]
 ***********************************************************************/

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
