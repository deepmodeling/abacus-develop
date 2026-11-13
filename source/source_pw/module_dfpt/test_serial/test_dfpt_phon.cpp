// Serial unit tests of the DFPT_Phon core (dfpt_phon.cpp): signed
// diagonalize against a closed-form spectrum, the isotropic LO-TO closed
// form, the Gamma check_sum_rule, the LO-TO direction normalization, the
// corrected-frequency diagonalization, and the output format regression.

#include "dfpt_test_tool.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_basis/module_pw/pw_basis_k.h"
#include "source_cell/qlist.h"
#include "source_cell/unitcell.h"
#include "source_pw/module_dfpt/dfpt_pw_data.h"
#include "source_base/complexmatrix.h"
#include "source_base/constants.h"
#include "source_cell/stru_fac.h"
#include "source_psi/psi.h"
#include "source_pw/module_dfpt/dfpt_kq_basis.h"
#include "source_pw/module_dfpt/dfpt_pert.h"
#include "source_pw/module_dfpt/dfpt_phon.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

// per-file self-contained fixture (gtest style): owns its objects and
// calls the dfpt_test free-function utilities
class DFPTPhonSerialTest : public ::testing::Test
{
  protected:
    ModuleBase::Matrix3 latvec_;
    UnitCell ucell_;
    ModulePW::PW_Basis pw_rho_;
    ModulePW::PW_Basis_K pw_wfc_;
    ModuleCell::QList qlist_;
    ModuleDFPT::DFPT_PW_Data data_;
    ModuleBase::Vector3<double> q_cart_;
    Structure_Factor sf_;
    ModuleDFPT::DFPT_Pert pert_;
    ModuleDFPT::DFPT_Phon phon_;

    // default (k, q) of this file: q is generic and k = -q so k+q = 0
    const ModuleBase::Vector3<double> k_d_{-0.13, 0.0, -0.07};
    const ModuleBase::Vector3<double> q_d_{0.13, 0.0, 0.07};

    void SetUp() override
    {
        dfpt_test::setup_cell(ucell_, latvec_);
        SetupPhon(k_d_, q_d_);
    }

    void TearDown() override
    {
        delete[] ucell_.atoms;
        ucell_.atoms = nullptr;
    }

    // (re)initialize the bases and the pert/phon wiring for a given (k, q)
    // pair; SetUp uses the default fixture values
    void SetupPhon(const ModuleBase::Vector3<double>& k_d, const ModuleBase::Vector3<double>& q_d)
    {
        dfpt_test::setup_bases(pw_rho_, pw_wfc_, qlist_, data_, ucell_, k_d, q_d, 2, q_cart_);
        pert_.init(ucell_, &pw_rho_, &pw_wfc_, sf_);
        phon_.init(ucell_, &pw_rho_, &pert_);
    }

    // independent Ry/bohr^2/amu -> cm^-1 conversion used by diagonalize
    double RyBohr2AmuToCm1() const
    {
        const double amu_kg = 1.66053906660e-27; // CODATA amu in kg
        return std::sqrt(ModuleBase::RYDBERG_SI / amu_kg) / (0.529177210903e-10 * 2.0 * ModuleBase::PI * 2.99792458e10);
    }

    // common setup of the isotropic loto closed-form tests: zero 6x6
    // dynamical matrix, eps_inf = 3I, Z*_1 = 1, Z*_2 = 2, and the
    // two-atom 12/4 mass table for the mass lookup
    void SetupIsotropicLoto()
    {
        data_.set_dynmat(0, ModuleBase::ComplexMatrix(6, 6, true));
        ModuleBase::matrix eps(3, 3, true);
        for (int d = 0; d < 3; ++d)
        {
            eps(d, d) = 3.0;
        }
        data_.set_dielectric(eps);
        ModuleBase::matrix z1(3, 3, true);
        ModuleBase::matrix z2(3, 3, true);
        z1(0, 0) = z1(1, 1) = z1(2, 2) = 1.0;
        z2(0, 0) = z2(1, 1) = z2(2, 2) = 2.0;
        data_.set_born(0, z1);
        data_.set_born(1, z2);
        dfpt_test::make_two_atom_cell(ucell_);
    }
};
#include <algorithm>
#include <cmath>
#include <complex>
#include <string>
#include <vector>

TEST_F(DFPTPhonSerialTest, DiagonalizeKnownMatrix)
{
    // 2-atom layout so nat3 = 6. Hermitian blocks with known closed-form
    // spectra: [[a, c], [conj(c), b]] has eigenvalues (a+b)/2
    // +- sqrt(((a-b)/2)^2 + |c|^2)
    dfpt_test::make_two_atom_cell(ucell_);
    const double lam[6] = {0.04, 0.09, 0.16, -0.02, 0.01, 0.1225}; // Ry/bohr^2/amu
    ModuleBase::ComplexMatrix dyn(6, 6, true);
    for (int i = 0; i < 6; ++i)
    {
        dyn(i, i) = std::complex<double>(lam[i], 0.0);
    }
    dyn(0, 1) = std::complex<double>(0.01, 0.02);
    dyn(1, 0) = std::conj(dyn(0, 1));
    dyn(2, 3) = std::complex<double>(-0.03, 0.005);
    dyn(3, 2) = std::conj(dyn(2, 3));
    data_.set_dynmat(0, dyn);
    phon_.diagonalize(0, data_);
    const std::vector<double> freq = data_.get_phon_freq(0);
    ASSERT_EQ(freq.size(), 6u);
    std::vector<double> expect;
    auto block = [&expect](double a, double b, std::complex<double> c) {
        const double mid = 0.5 * (a + b);
        const double rad = std::sqrt(std::pow(0.5 * (a - b), 2) + std::norm(c));
        expect.push_back(mid + rad);
        expect.push_back(mid - rad);
    };
    block(lam[0], lam[1], dyn(0, 1)); // coupled pair
    block(lam[2], lam[3], dyn(2, 3)); // coupled pair
    expect.push_back(lam[4]);         // untouched diagonal
    expect.push_back(lam[5]);
    for (double& e : expect)
    {
        const double s = (e >= 0.0) ? 1.0 : -1.0;
        e = s * std::sqrt(std::abs(e)) * RyBohr2AmuToCm1();
    }
    std::sort(expect.begin(), expect.end());
    std::vector<double> got = freq;
    std::sort(got.begin(), got.end());
    for (int i = 0; i < 6; ++i)
    {
        EXPECT_NEAR(got[i], expect[i], 1.0e-6 * std::abs(expect[i]));
    }
}

TEST_F(DFPTPhonSerialTest, AddLotoIsotropicClosedForm)
{
    // isotropic eps_inf = 3, Born charges Z*_1 = 1, Z*_2 = 2, masses 12/4
    SetupIsotropicLoto();

    const ModuleBase::Vector3<double> qhat(1.0, 0.0, 0.0);
    phon_.add_loto(qhat, data_);

    // closed form: D_NAC(0x,1x) = 4pi e2/Omega * 1*2/(3) / sqrt(12*4)
    const double expect = ModuleBase::FOUR_PI * ModuleBase::e2 / ucell_.omega / 3.0 * 2.0 / std::sqrt(48.0);
    const ModuleBase::ComplexMatrix dyn = data_.get_dynmat(0);
    EXPECT_NEAR(std::abs(dyn(0, 3) - std::complex<double>(expect, 0.0)), 0.0, 1.0e-12);
    EXPECT_NEAR(std::abs(dyn(3, 0) - std::complex<double>(expect, 0.0)), 0.0, 1.0e-12);
    // off-qhat elements untouched
    EXPECT_DOUBLE_EQ(std::abs(dyn(1, 4)), 0.0);
}

TEST_F(DFPTPhonSerialTest, CheckSumRuleAtGamma)
{
    // the sum rule is a Gamma-only statement: use a Gamma q point
    qlist_.kvec_d[0] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
    // zero dynamical matrix trivially satisfies the rule
    data_.set_dynmat(0, ModuleBase::ComplexMatrix(3, 3, true));
    EXPECT_TRUE(phon_.check_sum_rule(0, data_));
    // a uniform constant shift violates it
    ModuleBase::ComplexMatrix dyn(3, 3, true);
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            dyn(i, j) = std::complex<double>(0.1, 0.0);
        }
    }
    data_.set_dynmat(0, dyn);
    EXPECT_FALSE(phon_.check_sum_rule(0, data_));
}

TEST_F(DFPTPhonSerialTest, LotoDirNormalization)
{
    // default is the isotropic (1,1,1)/sqrt(3)
    const ModuleBase::Vector3<double> def = data_.get_loto_dir();
    const double inv = 1.0 / std::sqrt(3.0);
    EXPECT_NEAR(def.x, inv, 1.0e-12);
    EXPECT_NEAR(def.y, inv, 1.0e-12);
    EXPECT_NEAR(def.z, inv, 1.0e-12);
    // any non-null vector is normalized to a unit direction
    data_.set_loto_dir(ModuleBase::Vector3<double>(2.0, 0.0, 0.0));
    const ModuleBase::Vector3<double> x = data_.get_loto_dir();
    EXPECT_NEAR(x.x, 1.0, 1.0e-12);
    EXPECT_NEAR(x.y, 0.0, 1.0e-12);
    EXPECT_NEAR(x.z, 0.0, 1.0e-12);
    EXPECT_NEAR(std::sqrt(x * x), 1.0, 1.0e-12);
    // a null vector keeps the previous direction
    data_.set_loto_dir(ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    EXPECT_NEAR(data_.get_loto_dir().x, 1.0, 1.0e-12);
}

TEST_F(DFPTPhonSerialTest, DiagonalizeLotoClosedForm)
{
    // same isotropic fixture as AddLotoIsotropicClosedForm: zero dynamical
    // matrix + eps = 3I, Z*_1 = 1, Z*_2 = 2, masses 12/4, qhat = x.
    // add_loto fills BOTH the diagonal and cross xx elements with
    // pref = 4pi e2/Omega/3: (0x,0x) = pref/12, (3x,3x) = pref,
    // (0x,3x) = pref*2/sqrt(48); the 2x2 block
    // [[1/12, 2/sqrt48], [2/sqrt48, 1]]*pref has eigenvalues
    // {0, 13/12 * pref} (determinant 1/12 - 4/48 = 0), the yy/zz blocks
    // stay zero, so the spectrum is {13/12*pref, 0 x 5} in Ry/bohr^2/amu
    SetupIsotropicLoto();

    phon_.add_loto(ModuleBase::Vector3<double>(1.0, 0.0, 0.0), data_);
    phon_.diagonalize_loto(data_);

    const double pref = ModuleBase::FOUR_PI * ModuleBase::e2 / ucell_.omega / 3.0;
    const double expect = std::sqrt(13.0 / 12.0 * pref) * RyBohr2AmuToCm1();
    const std::vector<double> freq = data_.get_phon_freq_loto();
    ASSERT_EQ(freq.size(), static_cast<size_t>(6));
    // signed spectrum: one +expect, five zeros (sorted); the zeros carry
    // zheev roundoff of order sqrt(eps_mach * lambda_max) in frequency
    std::vector<double> sorted = freq;
    std::sort(sorted.begin(), sorted.end());
    EXPECT_NEAR(sorted.back(), expect, 1.0e-6 * std::abs(expect));
    for (int i = 0; i < 5; ++i)
    {
        EXPECT_NEAR(sorted[i], 0.0, 1.0e-5);
    }
}

TEST_F(DFPTPhonSerialTest, FormatReportsRegression)
{
    // fixture q = (0.13, 0, 0.07) direct; three crafted frequencies
    data_.set_phon_freq(0, std::vector<double>{-7.32457, 517.491, 0.0});
    const std::string qrep = phon_.format_q_report(0, data_);
    const std::string expect_q = " DFPT phonon frequencies at q #0 = (0.130000 0.000000 0.070000) "
                                 "(direct) in cm^-1:\n"
                                 "   mode   0 : -7.324570 cm^-1\n"
                                 "   mode   1 : 517.491000 cm^-1\n"
                                 "   mode   2 : 0.000000 cm^-1\n";
    EXPECT_EQ(qrep, expect_q);

    // LO-TO report: empty before the corrected frequencies exist
    EXPECT_TRUE(phon_.format_loto_report(data_).empty());
    data_.set_loto_dir(ModuleBase::Vector3<double>(0.0, 3.0, 0.0));
    data_.set_phon_freq_loto(std::vector<double>{0.0, 520.123456, 520.123457});
    const std::string lrep = phon_.format_loto_report(data_);
    const std::string expect_l = " DFPT LO-TO corrected frequencies at q #0 along q->0 direction "
                                 "(0.000000 1.000000 0.000000) in cm^-1:\n"
                                 "   mode   0 : 0.000000 cm^-1\n"
                                 "   mode   1 : 520.123456 cm^-1\n"
                                 "   mode   2 : 520.123457 cm^-1\n";
    EXPECT_EQ(lrep, expect_l);
}
