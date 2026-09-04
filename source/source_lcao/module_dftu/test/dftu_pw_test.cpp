#include "gtest/gtest.h"
#include <complex>
#include <vector>
#include "source_io/module_parameter/parameter.h"
#include "source_io/module_parameter/test_parameters.h"
#include "source_base/matrix.h"
#include "source_pw/module_pwdft/dftu_base_tools.h"

/***********************************************************************
 * Unit tests for DFT+U PW nspin=1/2/4 support (PR-2)
 *
 * Test targets:
 *   1. Energy weight logic: weight_eu and diag_coeff for nspin=1/2/4
 *   2. Becp index logic: different index formulas for nspin=1/2 vs nspin=4
 *   3. pot_onsite effective potential: cal_occ_pw pot_onsite calculation for all nspin modes
 *   4. Energy calculation: E_U accumulation with correct weights
 *   5. Occupation matrix accumulation from becp: the core loop of cal_occ_pw
 *   6. Multi-atom split layout: [all_up | all_dn] layout for nspin=2
 *   7. OnsitePsOp kernel: pot_onsite application to ps for npol=1
 *
 * Strategy: test energy weights and becp index logic as pure
 * arithmetic — no need to link against full ABACUS libraries.
 * set_occ_mat is tested via integration tests.
 ***********************************************************************/

class DftuPwTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}
};

// =====================================================================
// Energy weight + becp index tests (merged from 5 tests)
// =====================================================================

TEST_F(DftuPwTest, EnergyWeightsAllNspin)
{
    struct Case { int nspin; double expected_weight; double expected_diag; };
    Case cases[] = {{1, 1.0, 0.5}, {2, 0.5, 0.5}, {4, 0.25, 1.0}};
    for (const auto& c : cases) {
        TestParameters::input().nspin = c.nspin;
        double weight_eu = 1;
        switch (PARAM.inp.nspin) {
            case 1: weight_eu = 1.0; break;
            case 2: weight_eu = 0.5; break;
            case 4: weight_eu = 0.25; break;
            default: break;
        }
        const double diag_coeff = PARAM.inp.nspin == 4 ? 1.0 : 0.5;
        EXPECT_DOUBLE_EQ(weight_eu, c.expected_weight);
        EXPECT_DOUBLE_EQ(diag_coeff, c.expected_diag);
    }
}

TEST_F(DftuPwTest, BecpIndexNspin12vs4)
{
    const int nkb = 10, begin_ih = 3, m_begin = 4, m = 2, ib = 5;
    // nspin=1/2: index = ib*nkb + begin_ih + m_begin + m
    const int idx12 = ib * nkb + begin_ih + m_begin + m;
    EXPECT_EQ(idx12, 59);
    // nspin=4: index = ib*2*nkb + begin_ih + m_begin + m
    const int idx4 = ib * 2 * nkb + begin_ih + m_begin + m;
    EXPECT_EQ(idx4, 109);
    EXPECT_NE(idx12, idx4);
}

// =====================================================================
// pot_onsite effective potential tests (cal_occ_pw logic)
// =====================================================================

TEST_F(DftuPwTest, PotOnsitePotNspin1_DiagonalLocale)
{
    // For nspin=1: pot_onsite[m1,m2] = U * (0.5*delta(m1,m2) - occ_mat[m2*m_size+m1])
    // With diagonal occ_mat: occ_mat[m,m] = 0.3
    const double U_val = 4.0;
    const int m_size = 5; // d-orbital: 2*2+1
    const int size = m_size * m_size;

    std::vector<double> occ_mat_c(size, 0.0);
    for (int m = 0; m < m_size; m++)
        occ_mat_c[m * m_size + m] = 0.3; // diagonal

    std::vector<std::complex<double>> pot_onsite(size, {0.0, 0.0});
    DFTU_BASE::compute_pot_onsite_scalar(pot_onsite.data(), occ_mat_c.data(), U_val, 0.5, 1.0, m_size);

    // diagonal: U*(0.5 - 0.3) = 4.0*0.2 = 0.8
    for (int m = 0; m < m_size; m++)
        EXPECT_DOUBLE_EQ(pot_onsite[m * m_size + m].real(), 0.8);
    // off-diagonal: U*(0 - 0) = 0
    EXPECT_DOUBLE_EQ(pot_onsite[0 * m_size + 1].real(), 0.0);
    EXPECT_DOUBLE_EQ(pot_onsite[1 * m_size + 0].real(), 0.0);
}

TEST_F(DftuPwTest, PotOnsitePotNspin2_TwoSpinChannels)
{
    // nspin=2: two independent spin channels with same formula pot_onsite = U*(0.5*delta - occ_mat)
    const double U_val = 5.0;
    const int m_size = 3;
    const int size = m_size * m_size;

    std::vector<double> occ_mat_up(size, 0.0);
    std::vector<double> occ_mat_dn(size, 0.0);
    occ_mat_up[0] = 0.4; // occ_mat_up(0,0) = 0.4
    occ_mat_dn[0] = 0.1; // occ_mat_dn(0,0) = 0.1

    std::vector<std::complex<double>> pot_onsite_up(size, {0.0, 0.0});
    std::vector<std::complex<double>> pot_onsite_dn(size, {0.0, 0.0});
    DFTU_BASE::compute_pot_onsite_scalar(pot_onsite_up.data(), occ_mat_up.data(), U_val, 0.5, 0.5, m_size);
    DFTU_BASE::compute_pot_onsite_scalar(pot_onsite_dn.data(), occ_mat_dn.data(), U_val, 0.5, 0.5, m_size);

    // pot_onsite_up[0,0] = U*(0.5 - 0.4) = 0.5
    EXPECT_DOUBLE_EQ(pot_onsite_up[0].real(), 0.5);
    // pot_onsite_dn[0,0] = U*(0.5 - 0.1) = 2.0
    EXPECT_DOUBLE_EQ(pot_onsite_dn[0].real(), 2.0);
}

TEST_F(DftuPwTest, PotOnsitePotNspin4_PauliTransform)
{
    // nspin=4: after computing pot_onsite in Pauli basis, transform to spin basis
    // pot_onsite_spin[0] = 0.5*(pot_onsite_pauli[0] + pot_onsite_pauli[3])
    // pot_onsite_spin[3] = 0.5*(pot_onsite_pauli[0] - pot_onsite_pauli[3])
    // pot_onsite_spin[1] = 0.5*(pot_onsite_pauli[1] + i*pot_onsite_pauli[2])
    // pot_onsite_spin[2] = 0.5*(pot_onsite_pauli[1] - i*pot_onsite_pauli[2])
    const int m_size = 1;
    const int size = m_size * m_size;

    // For a single (m1,m2) pair, test the Pauli->spin transform (in-place)
    std::complex<double> pot_onsite[4];
    pot_onsite[0] = {1.0, 0.0}; // charge channel
    pot_onsite[1] = {0.5, 0.0}; // sigma_x
    pot_onsite[2] = {0.3, 0.0}; // sigma_y
    pot_onsite[3] = {0.2, 0.0}; // sigma_z

    DFTU_BASE::pauli_to_spin_basis(pot_onsite, m_size);

    EXPECT_DOUBLE_EQ(pot_onsite[0].real(), 0.6);  // 0.5*(1.0+0.2)
    EXPECT_DOUBLE_EQ(pot_onsite[0].imag(), 0.0);
    EXPECT_DOUBLE_EQ(pot_onsite[3].real(), 0.4);  // 0.5*(1.0-0.2)
    EXPECT_DOUBLE_EQ(pot_onsite[3].imag(), 0.0);
    EXPECT_DOUBLE_EQ(pot_onsite[1].real(), 0.25); // 0.5*0.5
    EXPECT_DOUBLE_EQ(pot_onsite[1].imag(), 0.15); // 0.5*0.3
    EXPECT_DOUBLE_EQ(pot_onsite[2].real(), 0.25); // 0.5*0.5
    EXPECT_DOUBLE_EQ(pot_onsite[2].imag(), -0.15);// -0.5*0.3
}

// =====================================================================
// Energy calculation tests
// =====================================================================

TEST_F(DftuPwTest, EnergyNspin12_DiagonalLocale)
{
    // E_U = sum_{m1,m2} U * weight_eu * occ_mat[m2,m1] * occ_mat[m1,m2]
    // nspin=1: weight_eu = 1.0, nspin=2: weight_eu = 0.5
    const double U_val = 4.0;
    const int m_size = 3;
    const int size = m_size * m_size;

    std::vector<double> occ_mat_c(size, 0.0);
    occ_mat_c[0 * m_size + 0] = 0.5;
    occ_mat_c[1 * m_size + 1] = 0.3;
    occ_mat_c[2 * m_size + 2] = 0.2;

    // nspin=1: E = U * 1.0 * (0.5^2 + 0.3^2 + 0.2^2) = 4 * 0.38 = 1.52
    std::vector<std::complex<double>> pot_onsite_nspin1(size, {0.0, 0.0});
    double energy_u = DFTU_BASE::compute_pot_onsite_scalar(
        pot_onsite_nspin1.data(), occ_mat_c.data(), U_val, 0.5, 1.0, m_size);
    EXPECT_DOUBLE_EQ(energy_u, 1.52);

    // nspin=2: two spin channels, weight_eu = 0.5
    std::vector<double> occ_mat_up(size, 0.0), occ_mat_dn(size, 0.0);
    occ_mat_up[0] = 0.4; occ_mat_dn[0] = 0.6;
    std::vector<std::complex<double>> pot_onsite_up(size, {0.0, 0.0});
    std::vector<std::complex<double>> pot_onsite_dn(size, {0.0, 0.0});
    energy_u = 0.0;
    energy_u += DFTU_BASE::compute_pot_onsite_scalar(
        pot_onsite_up.data(), occ_mat_up.data(), U_val, 0.5, 0.5, m_size);
    energy_u += DFTU_BASE::compute_pot_onsite_scalar(
        pot_onsite_dn.data(), occ_mat_dn.data(), U_val, 0.5, 0.5, m_size);
    // E = U*0.5*(0.4^2 + 0.6^2) = 4*0.5*(0.16+0.36) = 1.04
    EXPECT_DOUBLE_EQ(energy_u, 1.04);
}

TEST_F(DftuPwTest, EnergyNspin4_WithOffDiagonal)
{
    // nspin=4: weight_eu = 0.25, includes off-diagonal Pauli components
    const double U_val = 2.0;
    const int m_size = 2;
    const int size = m_size * m_size;
    const double weight_eu = 0.25;

    // 4 Pauli components stored contiguously
    std::vector<double> occ_mat_c(size * 4, 0.0);
    // charge channel (is=0)
    occ_mat_c[0] = 0.5; occ_mat_c[1] = 0.1;
    occ_mat_c[2] = 0.1; occ_mat_c[3] = 0.5;
    // sigma_x (is=1)
    occ_mat_c[size + 0] = 0.2; occ_mat_c[size + 1] = 0.0;
    occ_mat_c[size + 2] = 0.0; occ_mat_c[size + 3] = 0.2;

    std::vector<std::complex<double>> pot_onsite(size * 4, {0.0, 0.0});
    double energy_u = DFTU_BASE::compute_pot_onsite_spinor(
        pot_onsite.data(), occ_mat_c.data(), U_val, 1.0, weight_eu, m_size);

    // is=0: 2*0.25*(0.5*0.5 + 0.1*0.1 + 0.1*0.1 + 0.5*0.5) = 0.26
    // is=1: 2*0.25*(0.2*0.2 + 0 + 0 + 0.2*0.2) = 0.04
    // is=2,3: 0
    EXPECT_DOUBLE_EQ(energy_u, 0.30);
}

// =====================================================================
// Occupation matrix accumulation from becp (cal_occ_pw core loop)
// =====================================================================

TEST_F(DftuPwTest, LocaleAccumNspin12)
{
    // nspin=1/2: occ_mat[m1*m_size+m2] += weight * real(conj(becp[m1]) * becp[m2])
    const int m_size = 3, nkb = 5, begin_ih = 0, m_begin = 0, nbands = 2, ik = 0;

    std::vector<std::complex<double>> becp(nbands * nkb, {0.0, 0.0});
    becp[0 * nkb + 0] = {1.0, 0.0}; becp[0 * nkb + 1] = {0.0, 1.0}; becp[0 * nkb + 2] = {0.5, 0.5};
    becp[1 * nkb + 0] = {0.5, 0.0}; becp[1 * nkb + 1] = {0.5, -0.5}; becp[1 * nkb + 2] = {0.0, 1.0};

    ModuleBase::matrix wg(1, nbands);
    wg(0, 0) = 1.0;
    wg(0, 1) = 0.5;

    std::vector<double> occ_mat_c(m_size * m_size, 0.0);
    DFTU_BASE::accumulate_occ_scalar(
        occ_mat_c.data(), becp.data(), nbands, nkb,
        begin_ih, m_begin, m_size, wg, ik);

    // band0, w=1.0: occ_mat[0,0] = 1.0*|1|^2 = 1.0
    // band1, w=0.5: occ_mat[0,0] = 0.5*|0.5|^2 = 0.125
    EXPECT_DOUBLE_EQ(occ_mat_c[0], 1.125);

    // occ_mat[1,1]: band0 = 1.0*|i|^2 = 1.0, band1 = 0.5*|(0.5,-0.5)|^2 = 0.25
    EXPECT_DOUBLE_EQ(occ_mat_c[4], 1.25);
}

TEST_F(DftuPwTest, LocaleAccumNspin4_PauliComponents)
{
    // nspin=4: 4 Pauli components from becp with npol=2
    // occ[0] = w * conj(becp_up[m1]) * becp_up[m2]
    // occ[1] = w * conj(becp_up[m1]) * becp_dn[m2]
    // occ[2] = w * conj(becp_dn[m1]) * becp_up[m2]
    // occ[3] = w * conj(becp_dn[m1]) * becp_dn[m2]
    // occ_mat[ind] += (occ[0]+occ[3]).real()       -- charge
    // occ_mat[ind+size] += (occ[1]+occ[2]).real()   -- sigma_x
    // occ_mat[ind+2*size] += (occ[1]-occ[2]).imag() -- sigma_y
    // occ_mat[ind+3*size] += (occ[0]-occ[3]).real() -- sigma_z
    const int m_size = 1, nkb = 2, nbands = 1, npol = 2, ik = 0;

    std::vector<std::complex<double>> becp(nbands * 2 * nkb, {0.0, 0.0});
    becp[0] = {0.8, 0.0};       // becp_up[m=0]
    becp[nkb] = {0.0, 0.6};     // becp_dn[m=0]

    const int size = m_size * m_size;
    std::vector<double> occ_mat_c(size * 4, 0.0);

    ModuleBase::matrix wg(1, nbands);
    wg(0, 0) = 1.0;
    DFTU_BASE::accumulate_occ_spinor(
        occ_mat_c.data(), becp.data(), nbands, npol, nkb,
        0, 0, m_size, wg, ik);

    // becp_up = (0.8, 0), becp_dn = (0, 0.6)
    // occ[0] = 0.64, occ[1] = (0, 0.48), occ[2] = (0, -0.48), occ[3] = 0.36
    EXPECT_DOUBLE_EQ(occ_mat_c[0], 1.0);    // charge: (0.64+0.36).real = 1.0
    EXPECT_DOUBLE_EQ(occ_mat_c[1], 0.0);    // sigma_x: (occ1+occ2).real = 0
    EXPECT_DOUBLE_EQ(occ_mat_c[2], 0.96);   // sigma_y: (occ1-occ2).imag = 0.96
    EXPECT_DOUBLE_EQ(occ_mat_c[3], 0.28);   // sigma_z: (occ0-occ3).real = 0.28
}

// =====================================================================
// Multi-atom split layout test for nspin=2 (simplified P0-1 bug fix)
// Verifies that the split layout [all_up | all_dn] works correctly
// with multiple correlated atoms
// =====================================================================

TEST_F(DftuPwTest, MultiAtomSplitLayout_Nspin2)
{
    // 2 correlated atoms with d-orbital (l=2)
    const int nat = 2, m_size = 5, size = m_size * m_size;
    const int P = nat * size, total = P * 2, half_size = P;

    // pot_uterm_pw_index: split layout, each atom gets `size` entries
    std::vector<int> pot_uterm_pw_index = {0, size};

    // Simulate occ_mat values for both atoms
    std::vector<double> loc_up[2], loc_dn[2];
    for (int i = 0; i < 2; i++) {
        loc_up[i].assign(size, 0.0); loc_dn[i].assign(size, 0.0);
        for (int m = 0; m < m_size; m++) {
            loc_up[i][m * m_size + m] = 0.8 - i * 0.1;
            loc_dn[i][m * m_size + m] = 0.2 + i * 0.1;
        }
    }

    // --- Write to uom_array using split layout ---
    std::vector<double> uom_array(total, 0.0);
    for (int iat = 0; iat < nat; iat++)
        for (int mm = 0; mm < size; mm++) {
            uom_array[pot_uterm_pw_index[iat] + mm] = loc_up[iat][mm];
            uom_array[half_size + pot_uterm_pw_index[iat] + mm] = loc_dn[iat][mm];
        }

    // Verify split layout: first half = all spin-up, second half = all spin-down
    EXPECT_DOUBLE_EQ(uom_array[0], 0.8);                // atom 0 up diagonal
    EXPECT_DOUBLE_EQ(uom_array[size], 0.7);             // atom 1 up diagonal
    EXPECT_DOUBLE_EQ(uom_array[half_size], 0.2);        // atom 0 dn diagonal
    EXPECT_DOUBLE_EQ(uom_array[half_size + size], 0.3); // atom 1 dn diagonal

    // --- Read back and verify round-trip ---
    for (int iat = 0; iat < nat; iat++)
        EXPECT_DOUBLE_EQ(uom_array[pot_uterm_pw_index[iat]], loc_up[iat][0]);

    // --- pot_onsite values in split layout ---
    const double U_val = 5.0;
    const double diag_coeff = 0.5;
    std::vector<std::complex<double>> pot_uterm_pw(total, {0.0, 0.0});

    // atom 0 spin-up pot_onsite
    std::complex<double>* pot_onsite_up_0 = &pot_uterm_pw[0];
    pot_onsite_up_0[0] = U_val * (diag_coeff - loc_up[0][0]);
    // atom 0 spin-down pot_onsite (split layout: offset by half_size)
    std::complex<double>* pot_onsite_dn_0 = &pot_uterm_pw[half_size];
    pot_onsite_dn_0[0] = U_val * (diag_coeff - loc_dn[0][0]);

    EXPECT_DOUBLE_EQ(pot_onsite_up_0[0].real(), -1.5); // 5*(0.5-0.8)
    EXPECT_DOUBLE_EQ(pot_onsite_dn_0[0].real(), 1.5);  // 5*(0.5-0.2)

    // Verify no overlap between atoms in pot_onsite arrays
    std::complex<double>* pot_onsite_up_1 = &pot_uterm_pw[size];
    pot_onsite_up_1[0] = U_val * (diag_coeff - loc_up[1][0]);
    EXPECT_NE(pot_onsite_up_0[0], pot_onsite_up_1[0]);
}

// =====================================================================
// OnsitePsOp kernel test (simplified npol=1 branch)
// Tests the pot_onsite application to ps without full ABACUS integration
// =====================================================================

TEST_F(DftuPwTest, OnsitePsOpKernel_Nspin2_Npol1)
{
    // Simulate the npol=1 branch of onsite_ps_op kernel
    const int npm = 4, tnp = 10, orb_l = 2, tlp1 = 2 * orb_l + 1, nat = 2;

    // pot_onsite array: 2 atoms, each with tlp1*tlp1 = 25 elements
    std::vector<std::complex<double>> pot_onsite(nat * tlp1 * tlp1);
    for (size_t i = 0; i < pot_onsite.size(); i++)
        pot_onsite[i] = {static_cast<double>(i + 1), 0.0};

    // ip_m: maps each projector to m index within its atom
    std::vector<int> ip_m = {0, 1, 2, 3, 4, 0, 1, 2, 3, 4};
    std::vector<int> ip_iat = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1};
    std::vector<int> pot_onsite_begin_iat = {0, tlp1 * tlp1};

    // becp: npm * tnp
    std::vector<std::complex<double>> becp(npm * tnp, {0.0, 0.0});
    for (int ib = 0; ib < npm; ib++)
        for (int ip = 0; ip < tnp; ip++)
            becp[ib * tnp + ip] = {static_cast<double>(ib + ip + 1), 0.0};

    // ps: tnp * npm
    std::vector<std::complex<double>> ps(tnp * npm, {0.0, 0.0});

    // Kernel logic for npol=1 (EXACT copy from onsite_op.cpp)
    for (int ib = 0; ib < npm; ib++) {
        for (int ip = 0; ip < tnp; ip++) {
            int m1 = ip_m[ip];
            if (m1 < 0) continue;
            int iat = ip_iat[ip];
            const std::complex<double>* pot_onsite_iat = pot_onsite.data() + pot_onsite_begin_iat[iat];
            int ip2_begin = ip - m1, ip2_end = ip - m1 + tlp1;
            const int psind = ip * npm + ib;
            for (int ip2 = ip2_begin; ip2 < ip2_end; ip2++) {
                int m2 = ip_m[ip2];
                ps[psind] += pot_onsite_iat[m1 * tlp1 + m2] * becp[ib * tnp + ip2];
            }
        }
    }

    // Verify ps[0] (ib=0, ip=0, m1=0, iat=0)
    // ps[0] = sum_{ip2=0..4} pot_onsite[0*tlp1+ip_m[ip2]] * becp[0*tnp+ip2]
    std::complex<double> expected = {0.0, 0.0};
    for (int ip2 = 0; ip2 < tlp1; ip2++)
        expected += pot_onsite[ip2] * becp[ip2];
    EXPECT_DOUBLE_EQ(ps[0].real(), expected.real());
    EXPECT_DOUBLE_EQ(ps[0].imag(), expected.imag());
}
