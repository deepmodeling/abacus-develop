// Serial unit tests of the acoustic sum rule family DFPT_Phon::apply_asr
// (currently inside dfpt_phon.cpp; this file's name anticipates the
// Phase-2 source split that moves the ASR family into dfpt_phon_asr.cpp):
// simple/crystal enforcement, equivalent-atom redistribution, Born-charge
// sum rule, ASR+LO-TO interplay, and the directional-cancellation and
// cross-Cartesian regressions.

#include "dfpt_phon_serial_fixture.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>

namespace
{
/// per-atom PHYSICAL residual norm of a 2-atom (nat3=6) mass-reduced
/// dynamical matrix: R(a,alpha,beta) = sum_b sqrt(m_b) D(3a+alpha, 3b+beta)
/// at a FIXED direction beta, then combined -- opposite x/y residuals
/// cannot cancel
double row_sum_norm(const UnitCell& ucell, const ModuleBase::ComplexMatrix& dyn, int ia)
{
    double sqm[2] = {std::sqrt(ucell.atoms[ucell.iat2it[0]].mass),
                     std::sqrt(ucell.atoms[ucell.iat2it[1]].mass)};
    double s2 = 0.0;
    for (int da = 0; da < 3; ++da)
    {
        for (int db = 0; db < 3; ++db)
        {
            std::complex<double> s(0.0, 0.0);
            for (int ib = 0; ib < 2; ++ib)
            {
                s += sqm[ib] * dyn(3 * ia + da, 3 * ib + db);
            }
            s2 += std::norm(s);
        }
    }
    return std::sqrt(s2);
}

/// Hermitian 2-atom matrix with an asymmetric translational violation:
/// self blocks A1 != A2, off block B, plus an extra bump on atom 0
ModuleBase::ComplexMatrix make_violating_dyn()
{
    ModuleBase::ComplexMatrix dyn(6, 6, true);
    const double self1 = 0.05;
    const double self2 = 0.04;
    const double off = 0.02;
    const double bump = 0.012; // asymmetry: only on atom 0, x row
    for (int d = 0; d < 3; ++d)
    {
        dyn(d, d) = std::complex<double>(self1 + (d == 0 ? bump : 0.0), 0.0);
        dyn(3 + d, 3 + d) = std::complex<double>(self2, 0.0);
        dyn(d, 3 + d) = std::complex<double>(off, 0.01 * d);
        dyn(3 + d, d) = std::conj(dyn(d, 3 + d));
    }
    return dyn;
}
} // namespace

TEST_F(DFPTPhonSerialTest, AsrSimpleZeroesRowSums)
{
    MakeTwoAtomCell();
    // rewire the q list to Gamma (the fixture default q is not 0)
    SetupPhon(k_d_, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    const ModuleBase::ComplexMatrix dyn0 = make_violating_dyn();
    data_.set_dynmat(0, dyn0);
    phon_.apply_asr("simple", 0, data_);
    const ModuleBase::ComplexMatrix dyn = data_.get_dynmat(0);
    // legacy fix leaves only the tiny imaginary Hermitization residual
    const double scale = 0.06;
    EXPECT_LT(row_sum_norm(ucell_, dyn, 0) / scale, 1.0e-12);
    EXPECT_LT(row_sum_norm(ucell_, dyn, 1) / scale, 1.0e-12);
    // Hermiticity preserved by the re-symmetrization pass
    for (int i = 0; i < 6; ++i)
    {
        for (int j = i + 1; j < 6; ++j)
        {
            EXPECT_EQ(dyn(i, j), std::conj(dyn(j, i)));
        }
    }
}

TEST_F(DFPTPhonSerialTest, AsrCrystalZeroesRowSumsExactlyAndProjectsAcousticToZero)
{
    MakeTwoAtomCell();
    // rewire the q list to Gamma (the fixture default q is not 0)
    SetupPhon(k_d_, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    data_.set_dynmat(0, make_violating_dyn());
    phon_.apply_asr("crystal", 0, data_);
    const ModuleBase::ComplexMatrix dyn = data_.get_dynmat(0);
    EXPECT_LT(row_sum_norm(ucell_, dyn, 0), 1.0e-14);
    EXPECT_LT(row_sum_norm(ucell_, dyn, 1), 1.0e-14);
    for (int i = 0; i < 6; ++i)
    {
        for (int j = i + 1; j < 6; ++j)
        {
            EXPECT_EQ(dyn(i, j), std::conj(dyn(j, i)));
        }
    }
    // the three translation eigenvalues are projected out exactly: the
    // uniform-displacement subspace maps to zero under (1-P) D (1-P).
    // The matrix residual is ~1e-17; in frequency units the square root
    // rescaling lifts that to ~1e-5 cm^-1, still an exact zero physically
    phon_.diagonalize(0, data_);
    const std::vector<double> freq = data_.get_phon_freq(0);
    std::vector<double> byabs = freq;
    std::sort(byabs.begin(), byabs.end(), [](double a, double b) { return std::abs(a) < std::abs(b); });
    for (int m = 0; m < 3; ++m)
    {
        EXPECT_NEAR(byabs[m], 0.0, 1.0e-4);
    }
}

TEST_F(DFPTPhonSerialTest, AsrCrystalKeepsEquivalentAtomsEquivalent)
{
    MakeTwoAtomCell();
    // symmetry-equivalent atoms carry equal masses; make the two fixture
    // species mass-equal so the equivalence assertion is meaningful
    ucell_.atoms[1].mass = ucell_.atoms[0].mass;
    // rewire the q list to Gamma (the fixture default q is not 0)
    SetupPhon(k_d_, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    // symmetric input: identical self blocks, common off block, plus a
    // violation spread equally (row sums nonzero but atom-symmetric)
    ModuleBase::ComplexMatrix dyn(6, 6, true);
    const double self = 0.05;
    const double off = 0.02;
    for (int d = 0; d < 3; ++d)
    {
        dyn(d, d) = std::complex<double>(self, 0.0);
        dyn(3 + d, 3 + d) = std::complex<double>(self, 0.0);
        dyn(d, 3 + d) = std::complex<double>(off, 0.0);
        dyn(3 + d, d) = std::complex<double>(off, 0.0);
    }
    data_.set_dynmat(0, dyn);
    phon_.apply_asr("crystal", 0, data_);
    const ModuleBase::ComplexMatrix out = data_.get_dynmat(0);
    // crystal redistributes the mean uniformly: both self blocks stay equal
    // and both row sums vanish (simple would dump them on the diagonals,
    // also staying equal here -- the discriminating case is the asymmetric
    // one above, where only crystal shares the correction)
    for (int d = 0; d < 3; ++d)
    {
        EXPECT_EQ(out(d, d), out(3 + d, 3 + d));
    }
}

TEST_F(DFPTPhonSerialTest, AsrCrystalCorrectsBornChargeSumRule)
{
    MakeTwoAtomCell();
    // rewire the q list to Gamma (the fixture default q is not 0)
    SetupPhon(k_d_, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    data_.set_dynmat(0, make_violating_dyn());
    ModuleBase::matrix z1(3, 3, true);
    ModuleBase::matrix z2(3, 3, true);
    z1(0, 0) = 1.3;
    z1(1, 1) = 1.1;
    z1(2, 2) = 1.2;
    z2(0, 0) = -1.1;
    z2(1, 1) = -1.0;
    z2(2, 2) = -1.05; // sums: 0.2 / 0.1 / 0.15 != 0
    data_.set_born(0, z1);
    data_.set_born(1, z2);
    phon_.apply_asr("crystal", 0, data_);
    const ModuleBase::matrix o1 = data_.get_born(0);
    const ModuleBase::matrix o2 = data_.get_born(1);
    for (int a = 0; a < 3; ++a)
    {
        EXPECT_NEAR(o1(a, a) + o2(a, a), 0.0, 1.0e-14);
    }
    // shifted by the per-component mean: 0.1 / 0.05 / 0.075
    EXPECT_NEAR(o1(0, 0), 1.3 - 0.1, 1.0e-14);
    EXPECT_NEAR(o2(1, 1), -1.1 + 0.05, 1.0e-14);
}

TEST_F(DFPTPhonSerialTest, AsrCrystalPlusLotoKeepsAcousticExactlyZero)
{
    // regression of the mass-weighting bug: with UNequal masses the
    // translation constraints must live on the physical (sqrt(m)-weighted)
    // vectors, otherwise the LO-TO nonanalytic term -- which decouples from
    // the physical acoustic subspace only through sum_a Z*_a = 0 -- mixes
    // into the projected null space and lifts one acoustic frequency (the
    // bn_cubic Gamma run showed +36 cm^-1 before the fix)
    MakeTwoAtomCell(); // masses 12 / 4: unequal, the trigger
    SetupPhon(k_d_, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    data_.set_dynmat(0, make_violating_dyn());
    ModuleBase::matrix eps(3, 3, true);
    for (int d = 0; d < 3; ++d)
    {
        eps(d, d) = 3.0;
    }
    data_.set_dielectric(eps);
    ModuleBase::matrix z1(3, 3, true);
    ModuleBase::matrix z2(3, 3, true);
    z1(0, 0) = z1(1, 1) = z1(2, 2) = 1.3;
    z2(0, 0) = z2(1, 1) = z2(2, 2) = -1.1; // sums 0.2 / 0.2 / 0.2 != 0
    data_.set_born(0, z1);
    data_.set_born(1, z2);
    phon_.apply_asr("crystal", 0, data_);
    // LO-TO along (1,1,1)/sqrt(3): the acoustic triplet must stay at zero
    phon_.add_loto(ModuleBase::Vector3<double>(1.0, 1.0, 1.0), data_);
    phon_.diagonalize_loto(data_);
    const std::vector<double> freq = data_.get_phon_freq_loto();
    ASSERT_EQ(freq.size(), 6u);
    std::vector<double> byabs = freq;
    std::sort(byabs.begin(), byabs.end(), [](double a, double b) { return std::abs(a) < std::abs(b); });
    for (int m = 0; m < 3; ++m)
    {
        EXPECT_NEAR(byabs[m], 0.0, 1.0e-4);
    }
}

TEST_F(DFPTPhonSerialTest, AsrDirectionalCancellationDetected)
{
    // regression of the audit-metric fix: the block [[1,-1,0],[-1,1,0]]
    // has zero row sums (sum over ALL columns) but a NONZERO per-direction
    // residual D*t_x; the audit and apply_asr report must not miss it
    MakeTwoAtomCell();
    SetupPhon(k_d_, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    ModuleBase::ComplexMatrix dyn(6, 6, true);
    // atom0 self-block with the cancellation pattern
    dyn(0, 0) = std::complex<double>(1.0, 0.0);
    dyn(0, 1) = std::complex<double>(-1.0, 0.0);
    dyn(1, 0) = std::complex<double>(-1.0, 0.0);
    dyn(1, 1) = std::complex<double>(1.0, 0.0);
    // atom1 self-block mirrors it (keeps total charge-conservation flavor)
    dyn(3, 3) = std::complex<double>(1.0, 0.0);
    dyn(3, 4) = std::complex<double>(-1.0, 0.0);
    dyn(4, 3) = std::complex<double>(-1.0, 0.0);
    dyn(4, 4) = std::complex<double>(1.0, 0.0);
    data_.set_dynmat(0, dyn);

    // per-direction residual BEFORE any correction is nonzero
    const double residual_before = row_sum_norm(ucell_, dyn, 0);
    EXPECT_GT(residual_before, 0.1);

    // crystal ASR zeroes it
    ModuleBase::ComplexMatrix dyn_copy = dyn;
    data_.set_dynmat(0, dyn_copy);
    phon_.apply_asr("crystal", 0, data_);
    const ModuleBase::ComplexMatrix out = data_.get_dynmat(0);
    EXPECT_LT(row_sum_norm(ucell_, out, 0), 1.0e-14);
    EXPECT_LT(row_sum_norm(ucell_, out, 1), 1.0e-14);
}

TEST_F(DFPTPhonSerialTest, AsrSimpleCrossCartesianRegression)
{
    // regression of the simple-mode Hermitian fix: the cross-Cartesian
    // coupling D(0,4) = D(4,0) = delta (atom0-x <-> atom1-y) has no exact
    // Hermitian onsite-only fix; simple must still produce a final matrix
    // that is Hermitian, satisfies all three translation constraints
    // exactly, and has a real diagonal
    MakeTwoAtomCell(); // masses 12/4 (unequal)
    SetupPhon(k_d_, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    ModuleBase::ComplexMatrix dyn(6, 6, true);
    const double delta = 0.07;
    dyn(0, 4) = std::complex<double>(delta, 0.0);
    dyn(4, 0) = std::complex<double>(delta, 0.0);
    // also a diagonal self term so the matrix is not purely offsite
    dyn(0, 0) = std::complex<double>(0.05, 0.0);
    dyn(4, 4) = std::complex<double>(0.04, 0.0);
    dyn(2, 2) = std::complex<double>(0.03, 0.0);
    dyn(5, 5) = std::complex<double>(0.06, 0.0);
    data_.set_dynmat(0, dyn);
    phon_.apply_asr("simple", 0, data_);
    const ModuleBase::ComplexMatrix out = data_.get_dynmat(0);

    // translation constraints hold exactly (per-direction residuals)
    EXPECT_LT(row_sum_norm(ucell_, out, 0), 1.0e-14);
    EXPECT_LT(row_sum_norm(ucell_, out, 1), 1.0e-14);

    // Hermiticity
    for (int i = 0; i < 6; ++i)
    {
        for (int j = i + 1; j < 6; ++j)
        {
            EXPECT_EQ(out(i, j), std::conj(out(j, i))) << "i=" << i << " j=" << j;
        }
    }

    // diagonal is real
    for (int i = 0; i < 6; ++i)
    {
        EXPECT_EQ(out(i, i).imag(), 0.0) << "i=" << i;
    }
}
