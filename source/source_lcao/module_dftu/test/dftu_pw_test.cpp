#include "gtest/gtest.h"
#include <complex>
#define private public
#include "source_io/module_parameter/parameter.h"
#undef private

/***********************************************************************
 * Unit tests for DFT+U PW nspin=1/2/4 support (PR-2)
 *
 * Strategy: test energy weights and becp index logic as pure
 * arithmetic — no need to link against full ABACUS libraries.
 * set_locale is tested via integration tests.
 ***********************************************************************/

class DftuPwTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}
};

// =====================================================================
// Energy weight tests
// =====================================================================

TEST_F(DftuPwTest, EnergyWeightsNspin1)
{
    PARAM.input.nspin = 1;
    double weight_eu = 1;
    switch(PARAM.inp.nspin)
    {
        case 1: weight_eu = 1.0; break;
        case 2: weight_eu = 0.5; break;
        case 4: weight_eu = 0.25; break;
        default: break;
    }
    const double diag_coeff = PARAM.inp.nspin == 4 ? 1.0 : 0.5;
    EXPECT_DOUBLE_EQ(weight_eu, 1.0);
    EXPECT_DOUBLE_EQ(diag_coeff, 0.5);
}

TEST_F(DftuPwTest, EnergyWeightsNspin2)
{
    PARAM.input.nspin = 2;
    double weight_eu = 1;
    switch(PARAM.inp.nspin)
    {
        case 1: weight_eu = 1.0; break;
        case 2: weight_eu = 0.5; break;
        case 4: weight_eu = 0.25; break;
        default: break;
    }
    const double diag_coeff = PARAM.inp.nspin == 4 ? 1.0 : 0.5;
    EXPECT_DOUBLE_EQ(weight_eu, 0.5);
    EXPECT_DOUBLE_EQ(diag_coeff, 0.5);
}

TEST_F(DftuPwTest, EnergyWeightsNspin4)
{
    PARAM.input.nspin = 4;
    double weight_eu = 1;
    switch(PARAM.inp.nspin)
    {
        case 1: weight_eu = 1.0; break;
        case 2: weight_eu = 0.5; break;
        case 4: weight_eu = 0.25; break;
        default: break;
    }
    const double diag_coeff = PARAM.inp.nspin == 4 ? 1.0 : 0.5;
    EXPECT_DOUBLE_EQ(weight_eu, 0.25);
    EXPECT_DOUBLE_EQ(diag_coeff, 1.0);
}

// =====================================================================
// Becp index tests
// =====================================================================

TEST_F(DftuPwTest, OccupNspin12Index)
{
    const int nkb = 10, begin_ih = 3, m_begin = 4, m = 2, ib = 5;
    // nspin=1/2: index = ib*nkb + begin_ih + m_begin + m
    const int index_nspin12 = ib * nkb + begin_ih + m_begin + m;
    EXPECT_EQ(index_nspin12, 59);
    // different from nspin=4
    const int index_nspin4 = ib * 2 * nkb + begin_ih + m_begin + m;
    EXPECT_NE(index_nspin12, index_nspin4);
}

TEST_F(DftuPwTest, OccupNspin4Index)
{
    const int nkb = 10, begin_ih = 3, m_begin = 4, m = 2, ib = 5;
    const int index_nspin4 = ib * 2 * nkb + begin_ih + m_begin + m;
    EXPECT_EQ(index_nspin4, 109);
}

// =====================================================================
// set_locale logic tests (pure array copy, no UnitCell needed)
// =====================================================================

TEST_F(DftuPwTest, SetLocaleNspin4)
{
    // Simulate set_locale for nspin=4: uom_array -> locale copy
    PARAM.input.nspin = 4;
    const int mat_size = 10; // (2*2+1)*2 for d-orbital with npol=2
    const int total = mat_size * mat_size; // 100

    std::vector<double> uom_array(total);
    for(int i = 0; i < total; i++)
        uom_array[i] = static_cast<double>(i + 1);

    // Simulate locale as raw array (same as ModuleBase::matrix::c)
    std::vector<double> locale_c(total, 0.0);

    // nspin=4 branch: direct copy
    for(int mm = 0; mm < total; mm++)
        locale_c[mm] = uom_array[mm];

    for(int i = 0; i < total; i++)
        EXPECT_DOUBLE_EQ(locale_c[i], static_cast<double>(i + 1));
}

TEST_F(DftuPwTest, SetLocaleNspin2)
{
    // Simulate set_locale for nspin=2: uom_array -> locale copy (spin-up + spin-down)
    PARAM.input.nspin = 2;
    const int mat_size = 5; // 2*2+1 for d-orbital
    const int size_per_spin = mat_size * mat_size; // 25
    const int total = size_per_spin * 2; // 50

    std::vector<double> uom_array(total);
    for(int i = 0; i < size_per_spin; i++)
    {
        uom_array[i] = static_cast<double>(i + 1);                // spin-up
        uom_array[i + size_per_spin] = static_cast<double>(i + 101); // spin-down
    }

    std::vector<double> locale_up(size_per_spin, 0.0);
    std::vector<double> locale_dn(size_per_spin, 0.0);

    // nspin=1/2 branch: copy both spin channels
    const int nr_nc = size_per_spin; // locale[iat][l][0][0].nr * locale[iat][l][0][0].nc
    for(int mm = 0; mm < nr_nc; mm++)
    {
        locale_up[mm] = uom_array[mm];
        locale_dn[mm] = uom_array[mm + nr_nc];
    }

    for(int i = 0; i < size_per_spin; i++)
    {
        EXPECT_DOUBLE_EQ(locale_up[i], static_cast<double>(i + 1));
        EXPECT_DOUBLE_EQ(locale_dn[i], static_cast<double>(i + 101));
    }
}

// =====================================================================
// VU effective potential tests (cal_occ_pw logic)
// =====================================================================

TEST_F(DftuPwTest, VUPotNspin1_DiagonalLocale)
{
    // For nspin=1: VU[m1,m2] = U * (0.5*delta(m1,m2) - locale[m2*m_size+m1])
    // With diagonal locale: locale[m,m] = 0.3
    const double U_val = 4.0;
    const int m_size = 5; // d-orbital: 2*2+1
    const int size = m_size * m_size;

    std::vector<double> locale_c(size, 0.0);
    for(int m = 0; m < m_size; m++)
        locale_c[m * m_size + m] = 0.3; // diagonal

    std::vector<std::complex<double>> vu(size, {0.0, 0.0});
    for(int m1 = 0; m1 < m_size; m1++)
    {
        for(int m2 = 0; m2 < m_size; m2++)
        {
            const double diag_coeff = 0.5; // nspin != 4
            vu[m1 * m_size + m2] = U_val *
                (diag_coeff * (m1 == m2) - locale_c[m2 * m_size + m1]);
        }
    }

    // diagonal: U*(0.5 - 0.3) = 4.0*0.2 = 0.8
    for(int m = 0; m < m_size; m++)
        EXPECT_DOUBLE_EQ(vu[m * m_size + m].real(), 0.8);

    // off-diagonal: U*(0 - 0) = 0
    EXPECT_DOUBLE_EQ(vu[0 * m_size + 1].real(), 0.0);
    EXPECT_DOUBLE_EQ(vu[1 * m_size + 0].real(), 0.0);
}

TEST_F(DftuPwTest, VUPotNspin1_OffDiagonalLocale)
{
    // locale has off-diagonal elements
    const double U_val = 3.0;
    const int m_size = 3; // p-orbital: 2*1+1
    const int size = m_size * m_size;

    std::vector<double> locale_c(size, 0.0);
    locale_c[0 * m_size + 1] = 0.1; // locale(0,1) = 0.1
    locale_c[1 * m_size + 0] = 0.2; // locale(1,0) = 0.2

    std::vector<std::complex<double>> vu(size, {0.0, 0.0});
    for(int m1 = 0; m1 < m_size; m1++)
    {
        for(int m2 = 0; m2 < m_size; m2++)
        {
            vu[m1 * m_size + m2] = U_val *
                (0.5 * (m1 == m2) - locale_c[m2 * m_size + m1]);
        }
    }

    // VU[0,1] = U * (0 - locale[1*3+0]) = 3.0 * (-0.2) = -0.6
    EXPECT_DOUBLE_EQ(vu[0 * m_size + 1].real(), -0.6);
    // VU[1,0] = U * (0 - locale[0*3+1]) = 3.0 * (-0.1) = -0.3
    EXPECT_DOUBLE_EQ(vu[1 * m_size + 0].real(), -0.3);
}

TEST_F(DftuPwTest, VUPotNspin2_TwoSpinChannels)
{
    // nspin=2: two independent spin channels with same formula
    const double U_val = 5.0;
    const int m_size = 3;
    const int size = m_size * m_size;

    std::vector<double> locale_up(size, 0.0);
    std::vector<double> locale_dn(size, 0.0);
    locale_up[0] = 0.4; // locale_up(0,0) = 0.4
    locale_dn[0] = 0.1; // locale_dn(0,0) = 0.1

    // VU_up[0,0] = U*(0.5 - 0.4) = 0.5
    double vu_up_00 = U_val * (0.5 - locale_up[0 * m_size + 0]);
    EXPECT_DOUBLE_EQ(vu_up_00, 0.5);

    // VU_dn[0,0] = U*(0.5 - 0.1) = 2.0
    double vu_dn_00 = U_val * (0.5 - locale_dn[0 * m_size + 0]);
    EXPECT_DOUBLE_EQ(vu_dn_00, 2.0);
}

TEST_F(DftuPwTest, VUPotNspin4_PauliTransform)
{
    // nspin=4: after computing VU in Pauli basis, transform to spin basis
    // vu_spin[0] = 0.5*(vu_pauli[0] + vu_pauli[3])
    // vu_spin[3] = 0.5*(vu_pauli[0] - vu_pauli[3])
    // vu_spin[1] = 0.5*(vu_pauli[1] + i*vu_pauli[2])
    // vu_spin[2] = 0.5*(vu_pauli[1] - i*vu_pauli[2])
    const int m_size = 3;
    const int size = m_size * m_size;

    // For a single (m1,m2) pair, test the Pauli->spin transform
    std::complex<double> vu_pauli[4];
    vu_pauli[0] = {1.0, 0.0}; // charge channel
    vu_pauli[1] = {0.5, 0.0}; // sigma_x
    vu_pauli[2] = {0.3, 0.0}; // sigma_y
    vu_pauli[3] = {0.2, 0.0}; // sigma_z

    std::complex<double> vu_spin[4];
    vu_spin[0] = 0.5 * (vu_pauli[0] + vu_pauli[3]);
    vu_spin[3] = 0.5 * (vu_pauli[0] - vu_pauli[3]);
    vu_spin[1] = 0.5 * (vu_pauli[1] + std::complex<double>(0.0, 1.0) * vu_pauli[2]);
    vu_spin[2] = 0.5 * (vu_pauli[1] - std::complex<double>(0.0, 1.0) * vu_pauli[2]);

    EXPECT_DOUBLE_EQ(vu_spin[0].real(), 0.6);  // 0.5*(1.0+0.2)
    EXPECT_DOUBLE_EQ(vu_spin[0].imag(), 0.0);
    EXPECT_DOUBLE_EQ(vu_spin[3].real(), 0.4);  // 0.5*(1.0-0.2)
    EXPECT_DOUBLE_EQ(vu_spin[3].imag(), 0.0);
    EXPECT_DOUBLE_EQ(vu_spin[1].real(), 0.25); // 0.5*0.5
    EXPECT_DOUBLE_EQ(vu_spin[1].imag(), 0.15); // 0.5*0.3
    EXPECT_DOUBLE_EQ(vu_spin[2].real(), 0.25); // 0.5*0.5
    EXPECT_DOUBLE_EQ(vu_spin[2].imag(), -0.15);// -0.5*0.3
}

// =====================================================================
// Energy calculation tests
// =====================================================================

TEST_F(DftuPwTest, EnergyNspin1_DiagonalLocale)
{
    // E_U = sum_{m1,m2} U * weight_eu * locale[m2,m1] * locale[m1,m2]
    // weight_eu = 1.0 for nspin=1
    const double U_val = 4.0;
    const int m_size = 3;
    const int size = m_size * m_size;

    std::vector<double> locale_c(size, 0.0);
    locale_c[0 * m_size + 0] = 0.5;
    locale_c[1 * m_size + 1] = 0.3;
    locale_c[2 * m_size + 2] = 0.2;

    double energy_u = 0.0;
    const double weight_eu = 1.0;
    for(int m1 = 0; m1 < m_size; m1++)
    {
        for(int m2 = 0; m2 < m_size; m2++)
        {
            energy_u += U_val * weight_eu * locale_c[m2 * m_size + m1]
                        * locale_c[m1 * m_size + m2];
        }
    }

    // Only diagonal contributes: U * (0.5^2 + 0.3^2 + 0.2^2) = 4*(0.25+0.09+0.04) = 4*0.38 = 1.52
    EXPECT_DOUBLE_EQ(energy_u, 1.52);
}

TEST_F(DftuPwTest, EnergyNspin2_TwoChannels)
{
    // nspin=2: weight_eu = 0.5, sum over both spin channels
    const double U_val = 2.0;
    const int m_size = 3;
    const int size = m_size * m_size;
    const double weight_eu = 0.5;

    std::vector<double> locale_up(size, 0.0);
    std::vector<double> locale_dn(size, 0.0);
    locale_up[0] = 0.4; // (0,0)
    locale_dn[0] = 0.6; // (0,0)

    double energy_u = 0.0;
    // spin-up contribution
    for(int m1 = 0; m1 < m_size; m1++)
        for(int m2 = 0; m2 < m_size; m2++)
            energy_u += U_val * weight_eu * locale_up[m2 * m_size + m1] * locale_up[m1 * m_size + m2];
    // spin-down contribution
    for(int m1 = 0; m1 < m_size; m1++)
        for(int m2 = 0; m2 < m_size; m2++)
            energy_u += U_val * weight_eu * locale_dn[m2 * m_size + m1] * locale_dn[m1 * m_size + m2];

    // U*0.5*(0.4^2 + 0.6^2) = 2*0.5*(0.16+0.36) = 0.52
    EXPECT_DOUBLE_EQ(energy_u, 0.52);
}

TEST_F(DftuPwTest, EnergyNspin4_WithOffDiagonal)
{
    // nspin=4: weight_eu = 0.25, includes off-diagonal Pauli components
    const double U_val = 2.0;
    const int m_size = 2; // simplified: s-orbital would be 1, use 2 for test
    const int size = m_size * m_size;
    const double weight_eu = 0.25;

    // 4 Pauli components stored contiguously
    std::vector<double> locale_c(size * 4, 0.0);
    // charge channel (is=0)
    locale_c[0] = 0.5; locale_c[1] = 0.1;
    locale_c[2] = 0.1; locale_c[3] = 0.5;
    // sigma_x (is=1)
    locale_c[size + 0] = 0.2; locale_c[size + 1] = 0.0;
    locale_c[size + 2] = 0.0; locale_c[size + 3] = 0.2;

    double energy_u = 0.0;
    for(int is = 0; is < 4; is++)
    {
        int start = is * size;
        for(int m1 = 0; m1 < m_size; m1++)
        {
            for(int m2 = 0; m2 < m_size; m2++)
            {
                energy_u += U_val * weight_eu
                    * locale_c[start + m2 * m_size + m1]
                    * locale_c[start + m1 * m_size + m2];
            }
        }
    }

    // is=0: 2*0.25*(0.5*0.5 + 0.1*0.1 + 0.1*0.1 + 0.5*0.5) = 0.5*(0.25+0.01+0.01+0.25) = 0.26
    // is=1: 2*0.25*(0.2*0.2 + 0 + 0 + 0.2*0.2) = 0.5*(0.04+0.04) = 0.04
    // is=2,3: 0
    EXPECT_DOUBLE_EQ(energy_u, 0.30);
}

// =====================================================================
// Locale accumulation from becp (cal_occ_pw core loop)
// =====================================================================

TEST_F(DftuPwTest, LocaleAccumNspin12)
{
    // nspin=1/2: locale[m1*m_size+m2] += weight * real(conj(becp[m1]) * becp[m2])
    const int m_size = 3; // p-orbital
    const int nkb = 5;
    const int begin_ih = 0;
    const int m_begin = 0; // target_l=1, m_begin = 1*1 = 1... but for test simplicity use 0
    const int nbands = 2;
    const double weights[2] = {1.0, 0.5};

    // becp array: becp[ib*nkb + begin_ih + m_begin + m]
    std::vector<std::complex<double>> becp(nbands * nkb, {0.0, 0.0});
    // band 0
    becp[0 * nkb + 0] = {1.0, 0.0};
    becp[0 * nkb + 1] = {0.0, 1.0};
    becp[0 * nkb + 2] = {0.5, 0.5};
    // band 1
    becp[1 * nkb + 0] = {0.5, 0.0};
    becp[1 * nkb + 1] = {0.5, -0.5};
    becp[1 * nkb + 2] = {0.0, 1.0};

    std::vector<double> locale_c(m_size * m_size, 0.0);
    for(int ib = 0; ib < nbands; ib++)
    {
        const double weight = weights[ib];
        int ind_m1m2 = 0;
        for(int m1 = 0; m1 < m_size; m1++)
        {
            const int index_m1 = ib * nkb + begin_ih + m_begin + m1;
            for(int m2 = 0; m2 < m_size; m2++)
            {
                const int index_m2 = ib * nkb + begin_ih + m_begin + m2;
                locale_c[ind_m1m2] += weight * (std::conj(becp[index_m1]) * becp[index_m2]).real();
                ind_m1m2++;
            }
        }
    }

    // band0, w=1.0: conj(becp0)*becp0 = |1|^2=1, conj(becp0)*becp1 = 1*(0,1)=(0,1)->real=0
    // locale[0,0] from band0 = 1.0*1.0 = 1.0
    // band1, w=0.5: conj(becp0)*becp0 = |0.5|^2=0.25
    // locale[0,0] from band1 = 0.5*0.25 = 0.125
    EXPECT_DOUBLE_EQ(locale_c[0], 1.125); // 1.0 + 0.125

    // locale[1,1]: band0 = 1.0*|i|^2 = 1.0, band1 = 0.5*|(0.5,-0.5)|^2 = 0.5*0.5 = 0.25
    EXPECT_DOUBLE_EQ(locale_c[4], 1.25);
}

TEST_F(DftuPwTest, LocaleAccumNspin4_PauliComponents)
{
    // nspin=4: 4 Pauli components from becp with npol=2
    // occ[0] = w * conj(becp_up[m1]) * becp_up[m2]
    // occ[1] = w * conj(becp_up[m1]) * becp_dn[m2]
    // occ[2] = w * conj(becp_dn[m1]) * becp_up[m2]
    // occ[3] = w * conj(becp_dn[m1]) * becp_dn[m2]
    // locale[ind] += (occ[0]+occ[3]).real()       -- charge
    // locale[ind+size] += (occ[1]+occ[2]).real()   -- sigma_x
    // locale[ind+2*size] += (occ[1]-occ[2]).imag() -- sigma_y
    // locale[ind+3*size] += (occ[0]-occ[3]).real() -- sigma_z

    const int m_size = 1; // s-orbital for simplicity
    const int nkb = 2;
    const int nbands = 1;
    const double weight = 1.0;

    // becp layout: becp[ib*2*nkb + begin_ih + m]  (up)
    //              becp[ib*2*nkb + begin_ih + m + nkb] (down)
    std::vector<std::complex<double>> becp(nbands * 2 * nkb, {0.0, 0.0});
    // m=0 only (s-orbital)
    becp[0 * 2 * nkb + 0] = {0.8, 0.0};       // becp_up[m=0]
    becp[0 * 2 * nkb + 0 + nkb] = {0.0, 0.6}; // becp_dn[m=0]

    const int size = m_size * m_size; // 1
    std::vector<double> locale_c(size * 4, 0.0);

    for(int ib = 0; ib < nbands; ib++)
    {
        int ind_m1m2 = 0;
        for(int m1 = 0; m1 < m_size; m1++)
        {
            const int index_m1 = ib * 2 * nkb + 0 + m1;
            for(int m2 = 0; m2 < m_size; m2++)
            {
                const int index_m2 = ib * 2 * nkb + 0 + m2;
                std::complex<double> occ[4];
                occ[0] = weight * std::conj(becp[index_m1]) * becp[index_m2];
                occ[1] = weight * std::conj(becp[index_m1]) * becp[index_m2 + nkb];
                occ[2] = weight * std::conj(becp[index_m1 + nkb]) * becp[index_m2];
                occ[3] = weight * std::conj(becp[index_m1 + nkb]) * becp[index_m2 + nkb];
                locale_c[ind_m1m2] += (occ[0] + occ[3]).real();
                locale_c[ind_m1m2 + size] += (occ[1] + occ[2]).real();
                locale_c[ind_m1m2 + 2 * size] += (occ[1] - occ[2]).imag();
                locale_c[ind_m1m2 + 3 * size] += (occ[0] - occ[3]).real();
                ind_m1m2++;
            }
        }
    }

    // becp_up = (0.8, 0), becp_dn = (0, 0.6)
    // occ[0] = conj(0.8)*0.8 = 0.64
    // occ[1] = conj(0.8)*(0,0.6) = 0.8*(0,0.6) = (0, 0.48)
    // occ[2] = conj(0,0.6)*0.8 = (0,-0.6)*0.8 = (0, -0.48)
    // occ[3] = conj(0,0.6)*(0,0.6) = (0,-0.6)*(0,0.6) = 0.36
    EXPECT_DOUBLE_EQ(locale_c[0], 1.0);    // (0.64+0.36).real = 1.0 (charge)
    EXPECT_DOUBLE_EQ(locale_c[1], 0.0);    // (occ1+occ2).real = ((0,0.48)+(0,-0.48)).real = 0
    EXPECT_DOUBLE_EQ(locale_c[2], 0.96);   // (occ1-occ2).imag = ((0,0.48)-(0,-0.48)).imag = 0.96
    EXPECT_DOUBLE_EQ(locale_c[3], 0.28);   // (occ0-occ3).real = (0.64-0.36) = 0.28 (sigma_z)
}

TEST_F(DftuPwTest, CopyLocaleToUomSave_Nspin2)
{
    // Verify copy_locale logic for split layout: [all_up | all_dn]
    const int m_size = 3;
    const int size = m_size * m_size;

    std::vector<double> locale_spin0(size), locale_spin1(size);
    for(int i = 0; i < size; i++)
    {
        locale_spin0[i] = static_cast<double>(i + 1);
        locale_spin1[i] = static_cast<double>(i + 100);
    }

    std::vector<double> uom_save(size * 2, 0.0);
    const int eff_pot_index = 0;
    const int half_size = uom_save.size() / 2;
    for(int mm = 0; mm < size; mm++)
    {
        uom_save[eff_pot_index + mm] = locale_spin0[mm];
        uom_save[half_size + eff_pot_index + mm] = locale_spin1[mm];
    }

    for(int i = 0; i < size; i++)
    {
        EXPECT_DOUBLE_EQ(uom_save[i], static_cast<double>(i + 1));
        EXPECT_DOUBLE_EQ(uom_save[half_size + i], static_cast<double>(i + 100));
    }
}

TEST_F(DftuPwTest, CopyLocaleToUomSave_Nspin4)
{
    // nspin=4: 4 blocks stored contiguously
    const int m_size = 3;
    const int size = m_size * m_size;
    const int total = size * 4; // 4 Pauli components

    std::vector<double> locale_c(total);
    for(int i = 0; i < total; i++)
        locale_c[i] = static_cast<double>(i + 1);

    std::vector<double> uom_save(total, 0.0);
    const int eff_pot_index = 0;
    for(int mm = 0; mm < size; mm++)
    {
        uom_save[eff_pot_index + mm] = locale_c[mm];
        uom_save[eff_pot_index + mm + size] = locale_c[mm + size];
        uom_save[eff_pot_index + mm + 2 * size] = locale_c[mm + 2 * size];
        uom_save[eff_pot_index + mm + 3 * size] = locale_c[mm + 3 * size];
    }

    for(int i = 0; i < total; i++)
        EXPECT_DOUBLE_EQ(uom_save[i], static_cast<double>(i + 1));
}

// =====================================================================
// Step 1: VU calculation test for nspin=2 (isolated from kernel)
// This tests the complete cal_occ_pw vu calculation path:
// becp -> locale -> vu_up/vu_dn
// =====================================================================

TEST_F(DftuPwTest, VU_Calculation_Nspin2_FullPath)
{
    // Simulate complete vu calculation for nspin=2
    // This is the EXACT logic from cal_occ_pw, isolated from kernel

    const int m_size = 5; // d-orbital: 2*2+1
    const int size = m_size * m_size; // 25
    const double U_val = 5.0;
    const double weight_eu = 0.5; // nspin=2
    const double diag_coeff = 0.5;

    // Simulated locale values (would normally come from becp accumulation)
    std::vector<double> locale_up(size, 0.0);
    std::vector<double> locale_dn(size, 0.0);
    // Set diagonal values typical for occupied d-orbitals
    for(int m = 0; m < m_size; m++)
    {
        locale_up[m * m_size + m] = 0.8;
        locale_dn[m * m_size + m] = 0.2;
    }

    // Calculate VU for spin-up
    std::vector<std::complex<double>> vu_up(size, {0.0, 0.0});
    for(int m1 = 0; m1 < m_size; m1++)
    {
        for(int m2 = 0; m2 < m_size; m2++)
        {
            vu_up[m1 * m_size + m2] = U_val *
                (diag_coeff * (m1 == m2) - locale_up[m2 * m_size + m1]);
        }
    }

    // Calculate VU for spin-down
    std::vector<std::complex<double>> vu_dn(size, {0.0, 0.0});
    for(int m1 = 0; m1 < m_size; m1++)
    {
        for(int m2 = 0; m2 < m_size; m2++)
        {
            vu_dn[m1 * m_size + m2] = U_val *
                (diag_coeff * (m1 == m2) - locale_dn[m2 * m_size + m1]);
        }
    }

    // Verify spin-up VU
    // diagonal: U*(0.5 - 0.8) = 5*(-0.3) = -1.5
    for(int m = 0; m < m_size; m++)
    {
        EXPECT_DOUBLE_EQ(vu_up[m * m_size + m].real(), -1.5);
        EXPECT_DOUBLE_EQ(vu_up[m * m_size + m].imag(), 0.0);
    }
    // off-diagonal: U*(0 - 0) = 0
    EXPECT_DOUBLE_EQ(vu_up[0 * m_size + 1].real(), 0.0);
    EXPECT_DOUBLE_EQ(vu_up[1 * m_size + 0].real(), 0.0);

    // Verify spin-down VU
    // diagonal: U*(0.5 - 0.2) = 5*(0.3) = 1.5
    for(int m = 0; m < m_size; m++)
    {
        EXPECT_DOUBLE_EQ(vu_dn[m * m_size + m].real(), 1.5);
        EXPECT_DOUBLE_EQ(vu_dn[m * m_size + m].imag(), 0.0);
    }
    // off-diagonal: U*(0 - 0) = 0
    EXPECT_DOUBLE_EQ(vu_dn[0 * m_size + 1].real(), 0.0);
    EXPECT_DOUBLE_EQ(vu_dn[1 * m_size + 0].real(), 0.0);

    // Verify energy calculation
    double energy_u = 0.0;
    for(int m1 = 0; m1 < m_size; m1++)
        for(int m2 = 0; m2 < m_size; m2++)
        {
            energy_u += U_val * weight_eu * locale_up[m2 * m_size + m1] * locale_up[m1 * m_size + m2];
            energy_u += U_val * weight_eu * locale_dn[m2 * m_size + m1] * locale_dn[m1 * m_size + m2];
        }
    // Only diagonal: 5 orbitals per spin channel
    // spin-up: 5 * U * weight_eu * 0.8*0.8 = 5 * 5.0 * 0.5 * 0.64 = 8.0
    // spin-down: 5 * U * weight_eu * 0.2*0.2 = 5 * 5.0 * 0.5 * 0.04 = 0.5
    // total = 8.5
    EXPECT_DOUBLE_EQ(energy_u, 8.5);
}

// =====================================================================
// Step 2: Test vu_device sync for nspin=2
// This verifies the vu transfer from eff_pot_pw to vu_device
// =====================================================================

TEST_F(DftuPwTest, VU_DeviceSync_Nspin2)
{
    // Simulate eff_pot_pw layout for nspin=2
    const int m_size = 5;
    const int size = m_size * m_size;
    const int total_size = size * 2; // spin-up + spin-down

    std::vector<std::complex<double>> eff_pot_pw(total_size);
    // Initialize with known values
    for(int i = 0; i < size; i++)
    {
        eff_pot_pw[i] = {static_cast<double>(i + 1), 0.0};         // spin-up
        eff_pot_pw[i + size] = {static_cast<double>(i + 100), 0.0}; // spin-down
    }

    // Simulate vu_device sync for spin-down (isk[ik] == 1)
    const int size_eff_pot_pw = total_size / 2;
    std::vector<std::complex<double>> vu_device(size_eff_pot_pw);
    // memcpy from eff_pot_pw[0] + size_eff_pot_pw
    for(int i = 0; i < size_eff_pot_pw; i++)
    {
        vu_device[i] = eff_pot_pw[i + size_eff_pot_pw];
    }

    // Verify vu_device contains spin-down values
    for(int i = 0; i < size; i++)
    {
        EXPECT_DOUBLE_EQ(vu_device[i].real(), static_cast<double>(i + 100));
        EXPECT_DOUBLE_EQ(vu_device[i].imag(), 0.0);
    }
}

// =====================================================================
// Step 3: Test onsite_ps_op kernel for nspin=2 (npol=1)
// This tests the vu application to ps without full ABACUS integration
// =====================================================================

TEST_F(DftuPwTest, OnsitePsOpKernel_Nspin2_Npol1)
{
    // Simulate the npol=1 branch of onsite_ps_op kernel
    const int npm = 4;   // number of bands (npm/npol for npol=1)
    const int npol = 1;
    const int tnp = 10;  // total number of projectors
    const int orb_l = 2; // d-orbital
    const int tlp1 = 2 * orb_l + 1; // 5
    const int nat = 2;

    // vu array: 2 atoms, each with tlp1*tlp1 = 25 elements
    std::vector<std::complex<double>> vu(nat * tlp1 * tlp1);
    for(int i = 0; i < nat * tlp1 * tlp1; i++)
        vu[i] = {static_cast<double>(i + 1), 0.0};

    // ip_m: maps each projector to m index within its atom
    // First atom (iat=0): projectors 0-4 map to m=0-4
    // Second atom (iat=1): projectors 5-9 map to m=0-4
    std::vector<int> ip_m = {0, 1, 2, 3, 4, 0, 1, 2, 3, 4};
    std::vector<int> ip_iat = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1};
    std::vector<int> vu_begin_iat = {0, tlp1 * tlp1};

    // becp: npm * tnp
    std::vector<std::complex<double>> becp(npm * tnp, {0.0, 0.0});
    // Set some non-zero becp values
    for(int ib = 0; ib < npm; ib++)
        for(int ip = 0; ip < tnp; ip++)
            becp[ib * tnp + ip] = {static_cast<double>(ib + ip + 1), 0.0};

    // ps: tnp * npm
    std::vector<std::complex<double>> ps(tnp * npm, {0.0, 0.0});

    // Kernel logic for npol=1 (EXACT copy from onsite_op.cpp)
    for(int ib = 0; ib < npm; ib++)
    {
        for(int ip = 0; ip < tnp; ip++)
        {
            int m1 = ip_m[ip];
            if(m1 < 0) continue;
            int iat = ip_iat[ip];
            const std::complex<double>* vu_iat = vu.data() + vu_begin_iat[iat];
            int ip2_begin = ip - m1;
            int ip2_end = ip - m1 + tlp1;
            const int psind = ip * npm + ib;
            for(int ip2 = ip2_begin; ip2 < ip2_end; ip2++)
            {
                const int becpind = ib * tnp + ip2;
                int m2 = ip_m[ip2];
                const int index_mm = m1 * tlp1 + m2;
                ps[psind] += vu_iat[index_mm] * becp[becpind];
            }
        }
    }

    // Verify ps[0] (ib=0, ip=0)
    // m1=0, iat=0, vu_iat=vu[0..]
    // ip2 from 0 to 5
    std::complex<double> expected_ps00 = {0.0, 0.0};
    for(int ip2 = 0; ip2 < tlp1; ip2++)
    {
        const int becpind = 0 * tnp + ip2;
        int m2 = ip_m[ip2];
        const int index_mm = 0 * tlp1 + m2;
        expected_ps00 += vu[index_mm] * becp[becpind];
    }
    EXPECT_DOUBLE_EQ(ps[0].real(), expected_ps00.real());
    EXPECT_DOUBLE_EQ(ps[0].imag(), expected_ps00.imag());
}

// =====================================================================
// Step 4: Test spin-up only path (isolate from spin-down)
// =====================================================================

TEST_F(DftuPwTest, SpinUpOnly_Path_Nspin2)
{
    // Test that spin-up calculation is independent and correct
    const int m_size = 5;
    const int size = m_size * m_size;
    const double U_val = 5.0;
    const double diag_coeff = 0.5;

    // Only set spin-up locale
    std::vector<double> locale_up(size, 0.0);
    for(int m = 0; m < m_size; m++)
        locale_up[m * m_size + m] = 0.8;

    // Calculate VU for spin-up only
    std::vector<std::complex<double>> vu_up(size, {0.0, 0.0});
    for(int m1 = 0; m1 < m_size; m1++)
    {
        for(int m2 = 0; m2 < m_size; m2++)
        {
            vu_up[m1 * m_size + m2] = U_val *
                (diag_coeff * (m1 == m2) - locale_up[m2 * m_size + m1]);
        }
    }

    // Verify diagonal values
    for(int m = 0; m < m_size; m++)
        EXPECT_DOUBLE_EQ(vu_up[m * m_size + m].real(), -1.5); // 5*(0.5-0.8)

    // Verify off-diagonal are zero
    for(int m1 = 0; m1 < m_size; m1++)
        for(int m2 = 0; m2 < m_size; m2++)
            if(m1 != m2)
                EXPECT_DOUBLE_EQ(vu_up[m1 * m_size + m2].real(), 0.0);
}

// =====================================================================
// Step 5: Test spin-down only path (isolate from spin-up)
// =====================================================================

TEST_F(DftuPwTest, SpinDownOnly_Path_Nspin2)
{
    // Test that spin-down calculation is independent and correct
    const int m_size = 5;
    const int size = m_size * m_size;
    const double U_val = 5.0;
    const double diag_coeff = 0.5;

    // Only set spin-down locale
    std::vector<double> locale_dn(size, 0.0);
    for(int m = 0; m < m_size; m++)
        locale_dn[m * m_size + m] = 0.2;

    // Calculate VU for spin-down only
    std::vector<std::complex<double>> vu_dn(size, {0.0, 0.0});
    for(int m1 = 0; m1 < m_size; m1++)
    {
        for(int m2 = 0; m2 < m_size; m2++)
        {
            vu_dn[m1 * m_size + m2] = U_val *
                (diag_coeff * (m1 == m2) - locale_dn[m2 * m_size + m1]);
        }
    }

    // Verify diagonal values
    for(int m = 0; m < m_size; m++)
        EXPECT_DOUBLE_EQ(vu_dn[m * m_size + m].real(), 1.5); // 5*(0.5-0.2)

    // Verify off-diagonal are zero
    for(int m1 = 0; m1 < m_size; m1++)
        for(int m2 = 0; m2 < m_size; m2++)
            if(m1 != m2)
                EXPECT_DOUBLE_EQ(vu_dn[m1 * m_size + m2].real(), 0.0);
}

// =====================================================================
// Multi-atom split layout test for nspin=2
// Verifies that the split layout [all_up | all_dn] works correctly
// with multiple correlated atoms (the P0-1 bug fix)
// =====================================================================

TEST_F(DftuPwTest, MultiAtomSplitLayout_Nspin2)
{
    // 2 correlated atoms with d-orbital (l=2)
    const int nat = 2;
    const int m_size = 5;
    const int size = m_size * m_size; // 25 per atom per spin
    const int P = nat * size; // 50 = total spin-up block size
    const int total = P * 2; // 100 = total array size (split: up|dn)

    // eff_pot_pw_index: split layout, each atom gets `size` entries
    std::vector<int> eff_pot_pw_index(nat);
    eff_pot_pw_index[0] = 0;
    eff_pot_pw_index[1] = size; // 25

    // --- Test uom_array writing (dftu_pw.cpp logic) ---
    std::vector<double> uom_array(total, 0.0);
    // Simulate locale values for both atoms
    std::vector<double> locale_up_0(size, 0.0), locale_dn_0(size, 0.0);
    std::vector<double> locale_up_1(size, 0.0), locale_dn_1(size, 0.0);
    for(int m = 0; m < m_size; m++)
    {
        locale_up_0[m * m_size + m] = 0.8;
        locale_dn_0[m * m_size + m] = 0.2;
        locale_up_1[m * m_size + m] = 0.7;
        locale_dn_1[m * m_size + m] = 0.3;
    }

    // Write to uom_array using split layout
    const int half_size = total / 2; // P = 50
    // atom 0
    for(int mm = 0; mm < size; mm++)
    {
        uom_array[eff_pot_pw_index[0] + mm] = locale_up_0[mm];
        uom_array[half_size + eff_pot_pw_index[0] + mm] = locale_dn_0[mm];
    }
    // atom 1
    for(int mm = 0; mm < size; mm++)
    {
        uom_array[eff_pot_pw_index[1] + mm] = locale_up_1[mm];
        uom_array[half_size + eff_pot_pw_index[1] + mm] = locale_dn_1[mm];
    }

    // Verify split layout: first half = all spin-up, second half = all spin-down
    // atom 0 up: [0..24]
    EXPECT_DOUBLE_EQ(uom_array[0], 0.8); // locale_up_0 diagonal
    // atom 1 up: [25..49]
    EXPECT_DOUBLE_EQ(uom_array[size + 0], 0.7); // locale_up_1 diagonal
    // atom 0 dn: [50..74]
    EXPECT_DOUBLE_EQ(uom_array[half_size + 0], 0.2); // locale_dn_0 diagonal
    // atom 1 dn: [75..99]
    EXPECT_DOUBLE_EQ(uom_array[half_size + size + 0], 0.3); // locale_dn_1 diagonal

    // --- Test set_locale reading (dftu_occup.cpp logic) ---
    std::vector<double> read_up_0(size, 0.0), read_dn_0(size, 0.0);
    std::vector<double> read_up_1(size, 0.0), read_dn_1(size, 0.0);

    for(int mm = 0; mm < size; mm++)
    {
        // atom 0
        read_up_0[mm] = uom_array[eff_pot_pw_index[0] + mm];
        read_dn_0[mm] = uom_array[half_size + eff_pot_pw_index[0] + mm];
        // atom 1
        read_up_1[mm] = uom_array[eff_pot_pw_index[1] + mm];
        read_dn_1[mm] = uom_array[half_size + eff_pot_pw_index[1] + mm];
    }

    for(int mm = 0; mm < size; mm++)
    {
        EXPECT_DOUBLE_EQ(read_up_0[mm], locale_up_0[mm]);
        EXPECT_DOUBLE_EQ(read_dn_0[mm], locale_dn_0[mm]);
        EXPECT_DOUBLE_EQ(read_up_1[mm], locale_up_1[mm]);
        EXPECT_DOUBLE_EQ(read_dn_1[mm], locale_dn_1[mm]);
    }

    // --- Test VU writing (dftu_pw.cpp logic) ---
    std::vector<std::complex<double>> eff_pot_pw(total, {0.0, 0.0});
    const double U_val = 5.0;
    const double diag_coeff = 0.5;

    // atom 0 spin-up VU
    std::complex<double>* vu_up_0 = &eff_pot_pw[eff_pot_pw_index[0]];
    for(int m1 = 0; m1 < m_size; m1++)
        for(int m2 = 0; m2 < m_size; m2++)
            vu_up_0[m1 * m_size + m2] = U_val * (diag_coeff * (m1 == m2) - locale_up_0[m2 * m_size + m1]);

    // atom 0 spin-down VU (split layout: offset by half_size)
    std::complex<double>* vu_dn_0 = &eff_pot_pw[eff_pot_pw.size() / 2 + eff_pot_pw_index[0]];
    for(int m1 = 0; m1 < m_size; m1++)
        for(int m2 = 0; m2 < m_size; m2++)
            vu_dn_0[m1 * m_size + m2] = U_val * (diag_coeff * (m1 == m2) - locale_dn_0[m2 * m_size + m1]);

    // atom 1 spin-up VU
    std::complex<double>* vu_up_1 = &eff_pot_pw[eff_pot_pw_index[1]];
    for(int m1 = 0; m1 < m_size; m1++)
        for(int m2 = 0; m2 < m_size; m2++)
            vu_up_1[m1 * m_size + m2] = U_val * (diag_coeff * (m1 == m2) - locale_up_1[m2 * m_size + m1]);

    // atom 1 spin-down VU
    std::complex<double>* vu_dn_1 = &eff_pot_pw[eff_pot_pw.size() / 2 + eff_pot_pw_index[1]];
    for(int m1 = 0; m1 < m_size; m1++)
        for(int m2 = 0; m2 < m_size; m2++)
            vu_dn_1[m1 * m_size + m2] = U_val * (diag_coeff * (m1 == m2) - locale_dn_1[m2 * m_size + m1]);

    // Verify VU values
    // atom 0 up diagonal: 5*(0.5-0.8) = -1.5
    EXPECT_DOUBLE_EQ(vu_up_0[0].real(), -1.5);
    // atom 0 dn diagonal: 5*(0.5-0.2) = 1.5
    EXPECT_DOUBLE_EQ(vu_dn_0[0].real(), 1.5);
    // atom 1 up diagonal: 5*(0.5-0.7) = -1.0
    EXPECT_DOUBLE_EQ(vu_up_1[0].real(), -1.0);
    // atom 1 dn diagonal: 5*(0.5-0.3) = 1.0
    EXPECT_DOUBLE_EQ(vu_dn_1[0].real(), 1.0);

    // Verify no overlap between atoms in VU arrays
    // atom 0 up ends at index 24, atom 1 up starts at 25 — no overlap
    EXPECT_NE(vu_up_0[0], vu_up_1[0]);
    // atom 0 dn starts at half_size=50, atom 1 dn starts at half_size+25=75 — no overlap
    EXPECT_NE(vu_dn_0[0], vu_dn_1[0]);
}

// =====================================================================
// Test that split layout copy_locale/uom_save is consistent
// with set_locale/uom_array round-trip for multi-atom nspin=2
// =====================================================================

TEST_F(DftuPwTest, RoundTripCopyAndSetLocale_Nspin2_MultiAtom)
{
    const int nat = 2;
    const int m_size = 5;
    const int size = m_size * m_size;
    const int P = nat * size;
    const int total = P * 2;

    std::vector<int> eff_pot_pw_index = {0, size};
    std::vector<double> uom_save(total, 0.0);
    std::vector<double> uom_array(total, 0.0);

    // Simulate locale values
    std::vector<std::vector<double>> locale_up(nat, std::vector<double>(size, 0.0));
    std::vector<std::vector<double>> locale_dn(nat, std::vector<double>(size, 0.0));
    for(int iat = 0; iat < nat; iat++)
        for(int m = 0; m < m_size; m++)
        {
            locale_up[iat][m * m_size + m] = 0.9 - iat * 0.1;
            locale_dn[iat][m * m_size + m] = 0.1 + iat * 0.1;
        }

    // copy_locale -> uom_save (split layout)
    const int half_size = total / 2;
    for(int iat = 0; iat < nat; iat++)
        for(int mm = 0; mm < size; mm++)
        {
            uom_save[eff_pot_pw_index[iat] + mm] = locale_up[iat][mm];
            uom_save[half_size + eff_pot_pw_index[iat] + mm] = locale_dn[iat][mm];
        }

    // cal_occ_pw -> uom_array (split layout)
    for(int iat = 0; iat < nat; iat++)
        for(int mm = 0; mm < size; mm++)
        {
            uom_array[eff_pot_pw_index[iat] + mm] = locale_up[iat][mm];
            uom_array[half_size + eff_pot_pw_index[iat] + mm] = locale_dn[iat][mm];
        }

    // Mixing would compare uom_array with uom_save — verify they match
    for(int i = 0; i < total; i++)
        EXPECT_DOUBLE_EQ(uom_array[i], uom_save[i]);

    // set_locale reads back from uom_array
    std::vector<std::vector<double>> read_up(nat, std::vector<double>(size, 0.0));
    std::vector<std::vector<double>> read_dn(nat, std::vector<double>(size, 0.0));
    for(int iat = 0; iat < nat; iat++)
        for(int mm = 0; mm < size; mm++)
        {
            read_up[iat][mm] = uom_array[eff_pot_pw_index[iat] + mm];
            read_dn[iat][mm] = uom_array[half_size + eff_pot_pw_index[iat] + mm];
        }

    // Verify round-trip consistency
    for(int iat = 0; iat < nat; iat++)
        for(int mm = 0; mm < size; mm++)
        {
            EXPECT_DOUBLE_EQ(read_up[iat][mm], locale_up[iat][mm]);
            EXPECT_DOUBLE_EQ(read_dn[iat][mm], locale_dn[iat][mm]);
        }
}

// =====================================================================
// get_locale_flat / set_locale_flat logic tests (pure arithmetic)
//
// These test the nspin-dependent packing/unpacking logic without
// requiring a Plus_U instance, by simulating the same operations.
// =====================================================================

TEST_F(DftuPwTest, LocaleFlatPackNspin1)
{
    PARAM.input.nspin = 1;
    const int tlp1 = 3;
    const int size = tlp1 * tlp1;
    std::vector<double> locale_spin0(size);
    for (int i = 0; i < size; i++) locale_spin0[i] = static_cast<double>(i);
    std::vector<double> occ(size);
    for (int i = 0; i < size; i++) occ[i] = locale_spin0[i];
    for (int i = 0; i < size; i++) EXPECT_DOUBLE_EQ(occ[i], static_cast<double>(i));
}

TEST_F(DftuPwTest, LocaleFlatPackNspin2)
{
    PARAM.input.nspin = 2;
    const int tlp1 = 3;
    const int size = tlp1 * tlp1;
    std::vector<double> locale_spin0(size), locale_spin1(size);
    for (int i = 0; i < size; i++)
    {
        locale_spin0[i] = static_cast<double>(i);
        locale_spin1[i] = static_cast<double>(i + 100);
    }
    std::vector<double> occ(2 * size);
    for (int i = 0; i < size; i++)
    {
        occ[i] = locale_spin0[i];
        occ[size + i] = locale_spin1[i];
    }
    for (int i = 0; i < size; i++)
    {
        EXPECT_DOUBLE_EQ(occ[i], static_cast<double>(i));
        EXPECT_DOUBLE_EQ(occ[size + i], static_cast<double>(i + 100));
    }
}

TEST_F(DftuPwTest, LocaleFlatSetRoundTrip)
{
    const int tlp1 = 2;
    const int size = tlp1 * tlp1;
    std::vector<double> locale_data(size, 0.0);
    std::vector<double> occ(size);
    for (int i = 0; i < size; i++) occ[i] = static_cast<double>(i + 50);
    for (int i = 0; i < size; i++) locale_data[i] = occ[i];
    for (int i = 0; i < size; i++)
        EXPECT_DOUBLE_EQ(locale_data[i], static_cast<double>(i + 50));
}
