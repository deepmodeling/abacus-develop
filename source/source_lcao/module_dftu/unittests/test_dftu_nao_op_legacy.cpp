#include "gtest/gtest.h"
#include <cmath>
#include <complex>
#include <vector>

/***********************************************************************
 * Unit tests for DFT+U on-site potential and Pauli transfer
 * (dftu_nao_op_legacy.cpp).
 *
 * cal_pot_onsite: nspin=1,2: pot_onsite[is] = U*(0.5*delta - occ^T)
 *                 nspin=4:   pot_onsite[0] = U*(1.0*delta - occ^T),
 *                            pot_onsite[is>0] = -U*occ^T
 * transfer_pot_onsite: Pauli matrix transformation (nspin=4)
 ***********************************************************************/

// =====================================================================
// 1. cal_pot_onsite: Hubbard potential calculation
// =====================================================================

static void cal_pot_onsite(const std::vector<double>& occ, int m_size, double u_value,
                       std::vector<double>& pot_onsite, double& eu)
{
    pot_onsite.assign(occ.size(), 0.0);
    eu = 0.0;
    int spin_fold = occ.size() / m_size / m_size;
    if (spin_fold < 4) // nspin=1,2
    {
        for (int is = 0; is < spin_fold; ++is)
        {
            int start = is * m_size * m_size;
            for (int m1 = 0; m1 < m_size; m1++)
                for (int m2 = 0; m2 < m_size; m2++)
                {
                    pot_onsite[start + m1 * m_size + m2] = u_value * (0.5 * (m1 == m2) - occ[start + m2 * m_size + m1]);
                    eu += u_value * 0.5 * occ[start + m2 * m_size + m1] * occ[start + m1 * m_size + m2];
                }
        }
    }
    else // nspin=4
    {
        for (int m1 = 0; m1 < m_size; m1++)
            for (int m2 = 0; m2 < m_size; m2++)
            {
                pot_onsite[m1 * m_size + m2] = u_value * (1.0 * (m1 == m2) - occ[m2 * m_size + m1]);
                eu += u_value * 0.25 * occ[m2 * m_size + m1] * occ[m1 * m_size + m2];
            }
        for (int is = 1; is < spin_fold; ++is)
        {
            int start = is * m_size * m_size;
            for (int m1 = 0; m1 < m_size; m1++)
                for (int m2 = 0; m2 < m_size; m2++)
                {
                    pot_onsite[start + m1 * m_size + m2] = u_value * (0.0 - occ[start + m2 * m_size + m1]);
                    eu += u_value * 0.25 * occ[start + m2 * m_size + m1] * occ[start + m1 * m_size + m2];
                }
        }
    }
}

class CalVOfUTest : public ::testing::Test { protected: void SetUp() override {} };

TEST_F(CalVOfUTest, Nspin1_SingleOrbital_HalfFilled)
{
    std::vector<double> occ = {0.5};
    std::vector<double> pot_onsite; double eu = 0.0;
    cal_pot_onsite(occ, 1, 4.0, pot_onsite, eu);
    EXPECT_DOUBLE_EQ(pot_onsite[0], 0.0);
    EXPECT_DOUBLE_EQ(eu, 0.5);
}

TEST_F(CalVOfUTest, Nspin2_DOrbital_SpinPolarized)
{
    const int m_size = 5;
    std::vector<double> occ(m_size * m_size * 2, 0.0);
    for (int m = 0; m < m_size; m++) occ[m * m_size + m] = 0.8;
    for (int m = 0; m < m_size; m++) occ[m_size*m_size + m*m_size + m] = 0.2;
    std::vector<double> pot_onsite; double eu = 0.0;
    cal_pot_onsite(occ, m_size, 5.0, pot_onsite, eu);
    for (int m = 0; m < m_size; m++) EXPECT_NEAR(pot_onsite[m*m_size+m], -1.5, 1e-14);
    for (int m = 0; m < m_size; m++) EXPECT_NEAR(pot_onsite[m_size*m_size+m*m_size+m], 1.5, 1e-14);
    EXPECT_NEAR(eu, 8.5, 1e-14);
}

TEST_F(CalVOfUTest, Nspin4_Porbital_PauliBlocks)
{
    const int m_size = 3;
    std::vector<double> occ(m_size * m_size * 4, 0.0);
    for (int m = 0; m < m_size; m++) occ[m * m_size + m] = 0.5;
    std::vector<double> pot_onsite; double eu = 0.0;
    cal_pot_onsite(occ, m_size, 4.0, pot_onsite, eu);
    for (int m = 0; m < m_size; m++) EXPECT_NEAR(pot_onsite[m*m_size+m], 2.0, 1e-14);
    for (int is = 1; is < 4; is++)
        for (int i = 0; i < m_size*m_size; i++) EXPECT_NEAR(pot_onsite[is*m_size*m_size+i], 0.0, 1e-14);
    EXPECT_NEAR(eu, 0.75, 1e-14);
}

// =====================================================================
// 2. transfer_pot_onsite: Pauli matrix transformation (nspin=4)
// =====================================================================

static void transfer_pot_onsite(const std::vector<double>& pot_onsite_tmp,
                        std::vector<std::complex<double>>& pot_onsite)
{
    const int m_size = int(sqrt(pot_onsite_tmp.size()) / 2);
    const int m_size2 = m_size * m_size;
    pot_onsite.resize(pot_onsite_tmp.size());
    for (int m1 = 0; m1 < m_size; m1++)
        for (int m2 = 0; m2 < m_size; m2++)
        {
            int idx[4] = {m1*m_size+m2, m1*m_size+m2+m_size2, m2*m_size+m1+2*m_size2, m2*m_size+m1+3*m_size2};
            pot_onsite[idx[0]] = 0.5 * (pot_onsite_tmp[idx[0]] + pot_onsite_tmp[idx[3]]);
            pot_onsite[idx[3]] = 0.5 * (pot_onsite_tmp[idx[0]] - pot_onsite_tmp[idx[3]]);
            pot_onsite[idx[1]] = 0.5 * (pot_onsite_tmp[idx[1]] + std::complex<double>(0,1) * pot_onsite_tmp[idx[2]]);
            pot_onsite[idx[2]] = 0.5 * (pot_onsite_tmp[idx[1]] - std::complex<double>(0,1) * pot_onsite_tmp[idx[2]]);
        }
}

class Transferpot_onsiteTest : public ::testing::Test { protected: void SetUp() override {} };

TEST_F(Transferpot_onsiteTest, PauliI_IdentityInput)
{
    std::vector<double> pot_onsite_tmp = {1.0, 0.0, 0.0, 1.0};
    std::vector<std::complex<double>> pot_onsite;
    transfer_pot_onsite(pot_onsite_tmp, pot_onsite);
    EXPECT_NEAR(pot_onsite[0].real(), 1.0, 1e-15); EXPECT_NEAR(pot_onsite[0].imag(), 0.0, 1e-15);
    EXPECT_NEAR(pot_onsite[3].real(), 0.0, 1e-15);
}

TEST_F(Transferpot_onsiteTest, PureSigmaZ)
{
    std::vector<double> pot_onsite_tmp = {1.0, 0.0, 0.0, -1.0};
    std::vector<std::complex<double>> pot_onsite;
    transfer_pot_onsite(pot_onsite_tmp, pot_onsite);
    EXPECT_NEAR(pot_onsite[0].real(), 0.0, 1e-15); EXPECT_NEAR(pot_onsite[3].real(), 1.0, 1e-15);
}

TEST_F(Transferpot_onsiteTest, SigmaX_Y_Combined)
{
    std::vector<double> pot_onsite_tmp_x = {0.0, 1.0, 1.0, 0.0};
    std::vector<std::complex<double>> pot_onsite;
    transfer_pot_onsite(pot_onsite_tmp_x, pot_onsite);
    EXPECT_NEAR(pot_onsite[1].real(), 0.5, 1e-15); EXPECT_NEAR(pot_onsite[1].imag(), 0.5, 1e-15);
    EXPECT_NEAR(pot_onsite[2].real(), 0.5, 1e-15); EXPECT_NEAR(pot_onsite[2].imag(), -0.5, 1e-15);

    std::vector<double> pot_onsite_tmp_y = {0.0, 1.0, -1.0, 0.0};
    transfer_pot_onsite(pot_onsite_tmp_y, pot_onsite);
    EXPECT_NEAR(pot_onsite[1].real(), 0.5, 1e-15); EXPECT_NEAR(pot_onsite[1].imag(), -0.5, 1e-15);
    EXPECT_NEAR(pot_onsite[2].real(), 0.5, 1e-15); EXPECT_NEAR(pot_onsite[2].imag(), 0.5, 1e-15);
}
