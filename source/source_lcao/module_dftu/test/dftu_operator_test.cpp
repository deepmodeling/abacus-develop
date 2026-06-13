#include "gtest/gtest.h"
#include <cmath>
#include <complex>
#include <vector>
#include <unordered_map>

/***********************************************************************
 * Unit tests for DFT+U and DeltaSpin operator math and force/stress.
 * Tests cover: cal_v_of_u, transfer_vu, cal_coeff_lambda,
 * Force/IJR, Stress/IJR, Voigt->matrix, PW index setup
 ***********************************************************************/

// =====================================================================
// 1. cal_v_of_u: Hubbard potential calculation
// nspin=1,2 (spin_fold < 4): VU[is] = U * (0.5*delta - occ^T)
//   E_U += U * 0.5 * occ * occ^T
// nspin=4 (spin_fold == 4): VU[0] = U * (1.0*delta - occ^T), VU[is>0] = -U * occ^T
//   E_U += U * 0.25 * occ * occ^T
// =====================================================================

static void cal_v_of_u(const std::vector<double>& occ, int m_size, double u_value,
                       std::vector<double>& vu, double& eu)
{
    vu.assign(occ.size(), 0.0);
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
                    vu[start + m1 * m_size + m2] = u_value * (0.5 * (m1 == m2) - occ[start + m2 * m_size + m1]);
                    eu += u_value * 0.5 * occ[start + m2 * m_size + m1] * occ[start + m1 * m_size + m2];
                }
        }
    }
    else // nspin=4
    {
        // is=0: Pauli I block
        for (int m1 = 0; m1 < m_size; m1++)
            for (int m2 = 0; m2 < m_size; m2++)
            {
                vu[m1 * m_size + m2] = u_value * (1.0 * (m1 == m2) - occ[m2 * m_size + m1]);
                eu += u_value * 0.25 * occ[m2 * m_size + m1] * occ[m1 * m_size + m2];
            }
        // is=1,2,3: Pauli off-diagonal blocks
        for (int is = 1; is < spin_fold; ++is)
        {
            int start = is * m_size * m_size;
            for (int m1 = 0; m1 < m_size; m1++)
                for (int m2 = 0; m2 < m_size; m2++)
                {
                    vu[start + m1 * m_size + m2] = u_value * (0.0 - occ[start + m2 * m_size + m1]);
                    eu += u_value * 0.25 * occ[start + m2 * m_size + m1] * occ[start + m1 * m_size + m2];
                }
        }
    }
}

class CalVOfUTest : public ::testing::Test { protected: void SetUp() override {} };

TEST_F(CalVOfUTest, Nspin1_SingleOrbital_HalfFilled)
{
    std::vector<double> occ = {0.5};
    std::vector<double> vu; double eu = 0.0;
    cal_v_of_u(occ, 1, 4.0, vu, eu);
    // VU = 4*(0.5-0.5)=0, E_U = 4*0.5*0.5*0.5=0.5
    EXPECT_DOUBLE_EQ(vu[0], 0.0);
    EXPECT_DOUBLE_EQ(eu, 0.5);
}

TEST_F(CalVOfUTest, Nspin2_DOrbital_SpinPolarized)
{
    const int m_size = 5;
    std::vector<double> occ(m_size * m_size * 2, 0.0);
    for (int m = 0; m < m_size; m++) occ[m * m_size + m] = 0.8;        // spin-up majority
    for (int m = 0; m < m_size; m++) occ[m_size*m_size + m*m_size + m] = 0.2; // spin-down minority
    std::vector<double> vu; double eu = 0.0;
    cal_v_of_u(occ, m_size, 5.0, vu, eu);
    // spin-up VU: 5*(0.5-0.8)=-1.5, spin-down VU: 5*(0.5-0.2)=1.5
    for (int m = 0; m < m_size; m++) EXPECT_NEAR(vu[m*m_size+m], -1.5, 1e-14);
    for (int m = 0; m < m_size; m++) EXPECT_NEAR(vu[m_size*m_size+m*m_size+m], 1.5, 1e-14);
    // E_U = 5*0.5*[5*(0.8^2)+5*(0.2^2)] = 8.5
    EXPECT_NEAR(eu, 8.5, 1e-14);
}

TEST_F(CalVOfUTest, Nspin4_Porbital_PauliBlocks)
{
    const int m_size = 3;
    std::vector<double> occ(m_size * m_size * 4, 0.0);
    for (int m = 0; m < m_size; m++) occ[m * m_size + m] = 0.5; // Pauli I block
    std::vector<double> vu; double eu = 0.0;
    cal_v_of_u(occ, m_size, 4.0, vu, eu);
    // is=0: VU=4*(1.0-0.5)=2.0, is=1,2,3: VU=0
    for (int m = 0; m < m_size; m++) EXPECT_NEAR(vu[m*m_size+m], 2.0, 1e-14);
    for (int is = 1; is < 4; is++)
        for (int i = 0; i < m_size*m_size; i++) EXPECT_NEAR(vu[is*m_size*m_size+i], 0.0, 1e-14);
    // E_U = 4*0.25*3*(0.5*0.5) = 0.75
    EXPECT_NEAR(eu, 0.75, 1e-14);
}

// =====================================================================
// 2. transfer_vu: Pauli matrix transformation (nspin=4)
// vu[0] = 0.5*(vu_tmp[0]+vu_tmp[3])   // Pauli I
// vu[3] = 0.5*(vu_tmp[0]-vu_tmp[3])   // Pauli sigma_z
// vu[1] = 0.5*(vu_tmp[1]+i*vu_tmp[2]) // sigma_x+i*sigma_y
// vu[2] = 0.5*(vu_tmp[1]-i*vu_tmp[2]) // sigma_x-i*sigma_y
// =====================================================================

static void transfer_vu(const std::vector<double>& vu_tmp,
                        std::vector<std::complex<double>>& vu)
{
    const int m_size = int(sqrt(vu_tmp.size()) / 2);
    const int m_size2 = m_size * m_size;
    vu.resize(vu_tmp.size());
    for (int m1 = 0; m1 < m_size; m1++)
        for (int m2 = 0; m2 < m_size; m2++)
        {
            int idx[4] = {m1*m_size+m2, m1*m_size+m2+m_size2, m2*m_size+m1+2*m_size2, m2*m_size+m1+3*m_size2};
            vu[idx[0]] = 0.5 * (vu_tmp[idx[0]] + vu_tmp[idx[3]]);
            vu[idx[3]] = 0.5 * (vu_tmp[idx[0]] - vu_tmp[idx[3]]);
            vu[idx[1]] = 0.5 * (vu_tmp[idx[1]] + std::complex<double>(0,1) * vu_tmp[idx[2]]);
            vu[idx[2]] = 0.5 * (vu_tmp[idx[1]] - std::complex<double>(0,1) * vu_tmp[idx[2]]);
        }
}

class TransferVUTest : public ::testing::Test { protected: void SetUp() override {} };

TEST_F(TransferVUTest, PauliI_IdentityInput)
{
    std::vector<double> vu_tmp = {1.0, 0.0, 0.0, 1.0};
    std::vector<std::complex<double>> vu;
    transfer_vu(vu_tmp, vu);
    EXPECT_NEAR(vu[0].real(), 1.0, 1e-15); EXPECT_NEAR(vu[0].imag(), 0.0, 1e-15);
    EXPECT_NEAR(vu[3].real(), 0.0, 1e-15);
}

TEST_F(TransferVUTest, PureSigmaZ)
{
    std::vector<double> vu_tmp = {1.0, 0.0, 0.0, -1.0};
    std::vector<std::complex<double>> vu;
    transfer_vu(vu_tmp, vu);
    EXPECT_NEAR(vu[0].real(), 0.0, 1e-15); EXPECT_NEAR(vu[3].real(), 1.0, 1e-15);
}

TEST_F(TransferVUTest, SigmaX_Y_Combined)
{
    std::vector<double> vu_tmp_x = {0.0, 1.0, 1.0, 0.0};
    std::vector<std::complex<double>> vu;
    transfer_vu(vu_tmp_x, vu);
    EXPECT_NEAR(vu[1].real(), 0.5, 1e-15); EXPECT_NEAR(vu[1].imag(), 0.5, 1e-15);
    EXPECT_NEAR(vu[2].real(), 0.5, 1e-15); EXPECT_NEAR(vu[2].imag(), -0.5, 1e-15);

    std::vector<double> vu_tmp_y = {0.0, 1.0, -1.0, 0.0};
    transfer_vu(vu_tmp_y, vu);
    EXPECT_NEAR(vu[1].real(), 0.5, 1e-15); EXPECT_NEAR(vu[1].imag(), -0.5, 1e-15);
    EXPECT_NEAR(vu[2].real(), 0.5, 1e-15); EXPECT_NEAR(vu[2].imag(), 0.5, 1e-15);
}

// =====================================================================
// 3. cal_coeff_lambda: Lambda coefficient encoding
// Collinear (nspin=2): coeff[0]=lambda_z, coeff[1]=-lambda_z
// Non-collinear (nspin=4): coeff[0]=lambda_z, coeff[1]=lambda_x+i*lambda_y,
//   coeff[2]=lambda_x-i*lambda_y, coeff[3]=-lambda_z
// =====================================================================

static void cal_coeff_lambda_collinear(const std::vector<double>& lambda,
                                        std::vector<double>& coeff)
{ coeff[0] = lambda[0]; coeff[1] = -lambda[0]; }

static void cal_coeff_lambda_noncollinear(const std::vector<double>& lambda,
                                           std::vector<std::complex<double>>& coeff)
{
    coeff[0] = std::complex<double>(lambda[2], 0.0);
    coeff[1] = std::complex<double>(lambda[0], lambda[1]);
    coeff[2] = std::complex<double>(lambda[0], -lambda[1]);
    coeff[3] = std::complex<double>(-lambda[2], 0.0);
}

class CalCoeffLambdaTest : public ::testing::Test { protected: void SetUp() override {} };

TEST_F(CalCoeffLambdaTest, Collinear_PositiveLambdaZ)
{
    std::vector<double> lambda = {2.5}, coeff(2);
    cal_coeff_lambda_collinear(lambda, coeff);
    EXPECT_DOUBLE_EQ(coeff[0], 2.5); EXPECT_DOUBLE_EQ(coeff[1], -2.5);
}

TEST_F(CalCoeffLambdaTest, NonCollinear_General)
{
    std::vector<double> lambda = {1.0, 2.0, 3.0};
    std::vector<std::complex<double>> coeff(4);
    cal_coeff_lambda_noncollinear(lambda, coeff);
    EXPECT_NEAR(coeff[0].real(), 3.0, 1e-15);
    EXPECT_NEAR(coeff[1].real(), 1.0, 1e-15); EXPECT_NEAR(coeff[1].imag(), 2.0, 1e-15);
    EXPECT_NEAR(coeff[2].real(), 1.0, 1e-15); EXPECT_NEAR(coeff[2].imag(), -2.0, 1e-15);
    EXPECT_NEAR(coeff[3].real(), -3.0, 1e-15);
}

// =====================================================================
// 4. Force/IJR core loop
// force1 += VU * <d phi/dR1|chi> * <chi|phi> * DM
// force2 -= VU * <phi|chi> * <chi|phi> * DM
// nlm arrays: [value, deri_x, deri_y, deri_z]
// =====================================================================

static void cal_force_IJR_core(const std::vector<double>& vu_in,
    const std::vector<double>& nlm1, const std::vector<double>& nlm2,
    const double dm_val, int m_size, int nspin, double force1[3], double force2[3])
{
    const int m_size2 = m_size * m_size;
    for (int is = 0; is < nspin; is++)
        for (int m1 = 0; m1 < m_size; m1++)
            for (int m2 = 0; m2 < m_size; m2++)
            {
                double vu = vu_in[m1*m_size + m2 + is*m_size2], dm = dm_val;
                force1[0] += vu * nlm1[m1+m_size] * nlm2[m2] * dm;
                force1[1] += vu * nlm1[m1+2*m_size] * nlm2[m2] * dm;
                force1[2] += vu * nlm1[m1+3*m_size] * nlm2[m2] * dm;
                force2[0] -= vu * nlm1[m1+m_size] * nlm2[m2] * dm;
                force2[1] -= vu * nlm1[m1+2*m_size] * nlm2[m2] * dm;
                force2[2] -= vu * nlm1[m1+3*m_size] * nlm2[m2] * dm;
            }
}

class ForceIJRTest : public ::testing::Test { protected: void SetUp() override {} };

TEST_F(ForceIJRTest, SingleOrbital_SingleSpin)
{
    std::vector<double> vu = {2.0}, nlm1 = {1.0, 0.1, 0.2, 0.3}, nlm2 = {1.0, 0.0, 0.0, 0.0};
    double dm_val = 0.5, force1[3]={0}, force2[3]={0};
    cal_force_IJR_core(vu, nlm1, nlm2, dm_val, 1, 1, force1, force2);
    // force = VU*deri(nlm1)*val(nlm2)*DM = 2.0*{0.1,0.2,0.3}*1.0*0.5 = {0.1,0.2,0.3}
    EXPECT_NEAR(force1[0], 0.1, 1e-15); EXPECT_NEAR(force1[1], 0.2, 1e-15); EXPECT_NEAR(force1[2], 0.3, 1e-15);
    EXPECT_NEAR(force2[0], -0.1, 1e-15); EXPECT_NEAR(force2[1], -0.2, 1e-15); EXPECT_NEAR(force2[2], -0.3, 1e-15);
}

TEST_F(ForceIJRTest, ActionReaction)
{
    std::vector<double> vu = {1.5}, nlm1 = {1.0, 0.3, 0.4, 0.5}, nlm2 = {1.0, 0.0, 0.0, 0.0};
    double dm_val = 1.0, force1[3]={0}, force2[3]={0};
    cal_force_IJR_core(vu, nlm1, nlm2, dm_val, 1, 1, force1, force2);
    for (int i = 0; i < 3; i++) EXPECT_NEAR(force1[i], -force2[i], 1e-15);
}

// =====================================================================
// 5. Stress/IJR core loop
// stress[0]+=VU*DM*(nlm1_dx*dis1.x*nlm2_val+nlm1_val*nlm2_dx*dis2.x)
// stress[3]+=VU*DM*(nlm1_dy*dis1.y*nlm2_val+nlm1_val*nlm2_dy*dis2.y)
// stress[5]+=VU*DM*(nlm1_dz*dis1.z*nlm2_val+nlm1_val*nlm2_dz*dis2.z)
// =====================================================================

static void cal_stress_IJR_core(const std::vector<double>& vu_in,
    const std::vector<double>& nlm1, const std::vector<double>& nlm2,
    const double dm_val, int m_size, int nspin,
    double dis1[3], double dis2[3], double stress[6])
{
    const int m_size2 = m_size * m_size;
    for (int is = 0; is < nspin; is++)
        for (int m1 = 0; m1 < m_size; m1++)
            for (int m2 = 0; m2 < m_size; m2++)
            {
                double tmp = vu_in[m1*m_size+m2+is*m_size2] * dm_val;
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
    std::vector<double> vu = {1.0}, nlm1 = {1.0, 0.1, 0.0, 0.0}, nlm2 = {1.0, 0.2, 0.0, 0.0};
    double dm_val = 1.0, dis1[3] = {1.0, 0.0, 0.0}, dis2[3] = {-1.0, 0.0, 0.0}, stress[6] = {0.0};
    cal_stress_IJR_core(vu, nlm1, nlm2, dm_val, 1, 1, dis1, dis2, stress);
    // stress[0] = 1.0*(0.1*1.0*1.0 + 1.0*0.2*(-1.0)) = -0.1
    EXPECT_NEAR(stress[0], -0.1, 1e-15);
    EXPECT_NEAR(stress[1], 0.0, 1e-15); EXPECT_NEAR(stress[2], 0.0, 1e-15);
}

TEST_F(StressIJRTest, SymmetricDisplacement)
{
    std::vector<double> vu = {2.0}, nlm1 = {1.0, 0.1, 0.2, 0.3}, nlm2 = {1.0, 0.1, 0.2, 0.3};
    double dm_val = 1.0, dis1[3] = {1.0, 2.0, 3.0}, dis2[3] = {1.0, 2.0, 3.0}, stress[6] = {0.0};
    cal_stress_IJR_core(vu, nlm1, nlm2, dm_val, 1, 1, dis1, dis2, stress);
    EXPECT_NEAR(stress[0], 2.0*(0.1*1.0 + 1.0*0.1*1.0), 1e-15); // xx
    EXPECT_NEAR(stress[4], 2.0*(0.2*3.0 + 1.0*0.2*3.0), 1e-15); // yz
}

// =====================================================================
// 6. Stress Voigt -> matrix mapping
// Voigt [xx,xy,xz,yz,yy,zz] -> 3x3 symmetric matrix
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
    // [0,0]=1, [0,1]=2, [0,2]=3, [1,1]=5, [1,2]=4, [2,2]=6
    EXPECT_NEAR(matrix[0], 1.0, 1e-15); EXPECT_NEAR(matrix[1], 2.0, 1e-15);
    EXPECT_NEAR(matrix[2], 3.0, 1e-15); EXPECT_NEAR(matrix[3], 2.0, 1e-15);
    EXPECT_NEAR(matrix[4], 5.0, 1e-15); EXPECT_NEAR(matrix[5], 4.0, 1e-15);
    EXPECT_NEAR(matrix[6], 3.0, 1e-15); EXPECT_NEAR(matrix[7], 4.0, 1e-15);
    EXPECT_NEAR(matrix[8], 6.0, 1e-15);
    // Verify symmetry
    EXPECT_NEAR(matrix[1], matrix[3], 1e-15); EXPECT_NEAR(matrix[2], matrix[6], 1e-15);
    EXPECT_NEAR(matrix[5], matrix[7], 1e-15);
}

// =====================================================================
// 7. PW operator index setup (ip_iat, ip_m, vu_begin_iat)
// ip_m[ip] = m index if projector is correlated, else -1
// ip_iat[ip] = atom index, vu_begin_iat[iat] = VU array offset
// =====================================================================

class PWIndexSetupTest : public ::testing::Test
{
  protected:
    struct AtomInfo { int it, nh, target_l; }; // target_l=-1 if not correlated

    void setup_indices(const std::vector<AtomInfo>& atoms,
        std::vector<int>& ip_iat, std::vector<int>& ip_m,
        std::vector<int>& vu_begin_iat, int& vu_total_size)
    {
        int ip0 = 0, vu_begin = 0, npol = 1;
        ip_iat.resize(0); ip_m.resize(0); vu_begin_iat.resize(atoms.size());
        for (const auto& atom : atoms)
        {
            ip_iat.resize(ip_iat.size() + atom.nh);
            ip_m.resize(ip_m.size() + atom.nh);
            if (atom.target_l == -1)
            {
                for (int ip = 0; ip < atom.nh; ip++)
                { ip_iat[ip0] = static_cast<int>(&atom - &atoms[0]); ip_m[ip0++] = -1; }
                vu_begin_iat[&atom - &atoms[0]] = 0;
            }
            else
            {
                int tlp1 = 2 * atom.target_l + 1;
                vu_begin_iat[&atom - &atoms[0]] = vu_begin;
                vu_begin += tlp1 * tlp1 * npol * npol;
                int m_begin = atom.target_l * atom.target_l;
                int m_end = (atom.target_l + 1) * (atom.target_l + 1);
                for (int ip = 0; ip < atom.nh; ip++)
                {
                    ip_iat[ip0] = static_cast<int>(&atom - &atoms[0]);
                    ip_m[ip0++] = (ip >= m_begin && ip < m_end) ? ip - m_begin : -1;
                }
            }
        }
        vu_total_size = vu_begin;
    }
};

TEST_F(PWIndexSetupTest, SingleCorrelatedAtom_DOrbital)
{
    std::vector<AtomInfo> atoms = {{0, 9, 2}}; // s(1)+p(3)+d(5) projectors, l=2
    std::vector<int> ip_iat, ip_m, vu_begin_iat; int vu_total_size;
    setup_indices(atoms, ip_iat, ip_m, vu_begin_iat, vu_total_size);
    // Projectors 0-3 (s+p): m=-1; 4-8 (d): m=0..4
    EXPECT_EQ(ip_iat.size(), 9u);
    for (int ip = 0; ip < 4; ip++) EXPECT_EQ(ip_m[ip], -1);
    for (int ip = 4; ip < 9; ip++) { EXPECT_EQ(ip_iat[ip], 0); EXPECT_EQ(ip_m[ip], ip-4); }
    EXPECT_EQ(vu_begin_iat[0], 0);
    EXPECT_EQ(vu_total_size, 25); // 5*5
}

TEST_F(PWIndexSetupTest, MixedCorrelatedUncorrelated)
{
    std::vector<AtomInfo> atoms = {{0, 4, 1}, {1, 2, -1}}; // atom0: p-correlated, atom1: not
    std::vector<int> ip_iat, ip_m, vu_begin_iat; int vu_total_size;
    setup_indices(atoms, ip_iat, ip_m, vu_begin_iat, vu_total_size);
    // atom0: s(ip=0)->m=-1, p(ip=1,2,3)->m=0,1,2
    EXPECT_EQ(ip_iat[0], 0); EXPECT_EQ(ip_m[0], -1);
    EXPECT_EQ(ip_iat[1], 0); EXPECT_EQ(ip_m[1], 0);
    EXPECT_EQ(ip_iat[2], 0); EXPECT_EQ(ip_m[2], 1);
    EXPECT_EQ(ip_iat[3], 0); EXPECT_EQ(ip_m[3], 2);
    // atom1: all m=-1
    EXPECT_EQ(ip_iat[4], 1); EXPECT_EQ(ip_m[4], -1);
    EXPECT_EQ(ip_iat[5], 1); EXPECT_EQ(ip_m[5], -1);
    EXPECT_EQ(vu_total_size, 9); // 3*3 for p-orbital
}
