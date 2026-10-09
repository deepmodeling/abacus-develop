#include "gtest/gtest.h"
#include <cmath>
#include <vector>

/***********************************************************************
 * Unit tests for the force accumulation in dftu_nao_for_r.cpp.
 *
 * force1 += pot_onsite * <d phi/dR1|chi> * <chi|phi> * DM
 * force2 -= pot_onsite * <phi|chi> * <chi|phi> * DM
 * nlm arrays: [value, deri_x, deri_y, deri_z]
 ***********************************************************************/

static void cal_force_IJR_core(const std::vector<double>& pot_onsite_in,
    const std::vector<double>& nlm1, const std::vector<double>& nlm2,
    const double dm_val, int m_size, int nspin, double force1[3], double force2[3])
{
    const int m_size2 = m_size * m_size;
    for (int is = 0; is < nspin; is++)
        for (int m1 = 0; m1 < m_size; m1++)
            for (int m2 = 0; m2 < m_size; m2++)
            {
                double pot_onsite = pot_onsite_in[m1*m_size + m2 + is*m_size2], dm = dm_val;
                force1[0] += pot_onsite * nlm1[m1+m_size] * nlm2[m2] * dm;
                force1[1] += pot_onsite * nlm1[m1+2*m_size] * nlm2[m2] * dm;
                force1[2] += pot_onsite * nlm1[m1+3*m_size] * nlm2[m2] * dm;
                force2[0] -= pot_onsite * nlm1[m1+m_size] * nlm2[m2] * dm;
                force2[1] -= pot_onsite * nlm1[m1+2*m_size] * nlm2[m2] * dm;
                force2[2] -= pot_onsite * nlm1[m1+3*m_size] * nlm2[m2] * dm;
            }
}

class ForceIJRTest : public ::testing::Test { protected: void SetUp() override {} };

TEST_F(ForceIJRTest, SingleOrbital_SingleSpin)
{
    std::vector<double> pot_onsite = {2.0}, nlm1 = {1.0, 0.1, 0.2, 0.3}, nlm2 = {1.0, 0.0, 0.0, 0.0};
    double dm_val = 0.5, force1[3]={0}, force2[3]={0};
    cal_force_IJR_core(pot_onsite, nlm1, nlm2, dm_val, 1, 1, force1, force2);
    EXPECT_NEAR(force1[0], 0.1, 1e-15); EXPECT_NEAR(force1[1], 0.2, 1e-15); EXPECT_NEAR(force1[2], 0.3, 1e-15);
    EXPECT_NEAR(force2[0], -0.1, 1e-15); EXPECT_NEAR(force2[1], -0.2, 1e-15); EXPECT_NEAR(force2[2], -0.3, 1e-15);
}

TEST_F(ForceIJRTest, ActionReaction)
{
    std::vector<double> pot_onsite = {1.5}, nlm1 = {1.0, 0.3, 0.4, 0.5}, nlm2 = {1.0, 0.0, 0.0, 0.0};
    double dm_val = 1.0, force1[3]={0}, force2[3]={0};
    cal_force_IJR_core(pot_onsite, nlm1, nlm2, dm_val, 1, 1, force1, force2);
    for (int i = 0; i < 3; i++) EXPECT_NEAR(force1[i], -force2[i], 1e-15);
}
