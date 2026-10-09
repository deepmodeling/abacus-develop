#include "gtest/gtest.h"
#include <cmath>
#include <vector>

/***********************************************************************
 * Unit tests for the stress accumulation in dftu_nao_str_r.cpp.
 *
 * stress[0]+=pot_onsite*DM*(nlm1_dx*dis1.x*nlm2_val+nlm1_val*nlm2_dx*dis2.x)
 * stress[3]+=pot_onsite*DM*(nlm1_dy*dis1.y*nlm2_val+nlm1_val*nlm2_dy*dis2.y)
 * stress[5]+=pot_onsite*DM*(nlm1_dz*dis1.z*nlm2_val+nlm1_val*nlm2_dz*dis2.z)
 *
 * Also covers the Voigt -> 3x3 matrix mapping used by the stress output.
 ***********************************************************************/

// =====================================================================
// 1. Stress/IJR core loop
// =====================================================================

static void cal_stress_IJR_core(const std::vector<double>& pot_onsite_in,
    const std::vector<double>& nlm1, const std::vector<double>& nlm2,
    const double dm_val, int m_size, int nspin,
    double dis1[3], double dis2[3], double stress[6])
{
    const int m_size2 = m_size * m_size;
    for (int is = 0; is < nspin; is++)
        for (int m1 = 0; m1 < m_size; m1++)
            for (int m2 = 0; m2 < m_size; m2++)
            {
                double tmp = pot_onsite_in[m1*m_size+m2+is*m_size2] * dm_val;
                stress[0] += tmp*(nlm1[m1+m_size]*dis1[0]*nlm2[m2] + nlm1[m1]*nlm2[m2+m_size]*dis2[0]);
                stress[1] += tmp*(nlm1[m1+m_size]*dis1[1]*nlm2[m2] + nlm1[m1]*nlm2[m2+m_size]*dis2[1]);
                stress[2] += tmp*(nlm1[m1+m_size]*dis1[2]*nlm2[m2] + nlm1[m1]*nlm2[m2+m_size]*dis2[2]);
                stress[3] += tmp*(nlm1[m1+2*m_size]*dis1[1]*nlm2[m2] + nlm1[m1]*nlm2[m2+2*m_size]*dis2[1]);
                stress[4] += tmp*(nlm1[m1+2*m_size]*dis1[2]*nlm2[m2] + nlm1[m1]*nlm2[m2+2*m_size]*dis2[2]);
                stress[5] += tmp*(nlm1[m1+3*m_size]*dis1[2]*nlm2[m2] + nlm1[m1]*nlm2[m2+3*m_size]*dis2[2]);
            }
}

class StressIJRTest : public ::testing::Test { protected: void SetUp() override {} };

TEST_F(StressIJRTest, SingleOrbital_XDisplacement)
{
    std::vector<double> pot_onsite = {1.0}, nlm1 = {1.0, 0.1, 0.0, 0.0}, nlm2 = {1.0, 0.2, 0.0, 0.0};
    double dm_val = 1.0, dis1[3] = {1.0, 0.0, 0.0}, dis2[3] = {-1.0, 0.0, 0.0}, stress[6] = {0.0};
    cal_stress_IJR_core(pot_onsite, nlm1, nlm2, dm_val, 1, 1, dis1, dis2, stress);
    EXPECT_NEAR(stress[0], -0.1, 1e-15);
    EXPECT_NEAR(stress[1], 0.0, 1e-15); EXPECT_NEAR(stress[2], 0.0, 1e-15);
}

TEST_F(StressIJRTest, SymmetricDisplacement)
{
    std::vector<double> pot_onsite = {2.0}, nlm1 = {1.0, 0.1, 0.2, 0.3}, nlm2 = {1.0, 0.1, 0.2, 0.3};
    double dm_val = 1.0, dis1[3] = {1.0, 2.0, 3.0}, dis2[3] = {1.0, 2.0, 3.0}, stress[6] = {0.0};
    cal_stress_IJR_core(pot_onsite, nlm1, nlm2, dm_val, 1, 1, dis1, dis2, stress);
    EXPECT_NEAR(stress[0], 2.0*(0.1*1.0 + 1.0*0.1*1.0), 1e-15); // xx
    EXPECT_NEAR(stress[4], 2.0*(0.2*3.0 + 1.0*0.2*3.0), 1e-15); // yz
}

// =====================================================================
// 2. Stress Voigt -> matrix mapping
// =====================================================================

static void voigt_to_matrix(double stress_6[6], double matrix[9])
{
    for (int i = 0; i < 9; i++) matrix[i] = 0.0;
    matrix[0]=stress_6[0]; matrix[1]=stress_6[1]; matrix[2]=stress_6[2];
    matrix[3]=stress_6[1]; matrix[4]=stress_6[4]; matrix[5]=stress_6[3];
    matrix[6]=stress_6[2]; matrix[7]=stress_6[3]; matrix[8]=stress_6[5];
}

class VoigtToMatrixTest : public ::testing::Test { protected: void SetUp() override {} };

TEST_F(VoigtToMatrixTest, FullMappingAndSymmetry)
{
    double stress_6[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0}, matrix[9];
    voigt_to_matrix(stress_6, matrix);
    EXPECT_NEAR(matrix[0], 1.0, 1e-15); EXPECT_NEAR(matrix[1], 2.0, 1e-15);
    EXPECT_NEAR(matrix[2], 3.0, 1e-15); EXPECT_NEAR(matrix[3], 2.0, 1e-15);
    EXPECT_NEAR(matrix[4], 5.0, 1e-15); EXPECT_NEAR(matrix[5], 4.0, 1e-15);
    EXPECT_NEAR(matrix[6], 3.0, 1e-15); EXPECT_NEAR(matrix[7], 4.0, 1e-15);
    EXPECT_NEAR(matrix[8], 6.0, 1e-15);
    EXPECT_NEAR(matrix[1], matrix[3], 1e-15); EXPECT_NEAR(matrix[2], matrix[6], 1e-15);
    EXPECT_NEAR(matrix[5], matrix[7], 1e-15);
}
