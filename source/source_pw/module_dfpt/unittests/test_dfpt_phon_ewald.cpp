// Serial unit tests of DFPT_Phon::ion_ion (dfpt_phon_ewald.cpp):
// the Ewald force constants satisfy the acoustic sum rule at Gamma to
// grid accuracy (this pins the sign and magnitude of the Gaussian self
// term), give three zero acoustic modes at Gamma, and are cross-checked
// against a direct (unscreened) lattice Hessian sum at a generic
// incommensurate q where the dipole sum is oscillation-screened.

#include "dfpt_test_tool.h"
#include "source_base/complexmatrix.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_basis/module_pw/pw_basis_k.h"
#include "source_cell/qlist.h"
#include "source_cell/stru_fac.h"
#include "source_cell/unitcell.h"
#include "source_psi/psi.h"
#include "source_pw/module_dfpt/dfpt_kq_basis.h"
#include "source_pw/module_dfpt/dfpt_pert.h"
#include "source_pw/module_dfpt/dfpt_phon.h"
#include "source_pw/module_dfpt/dfpt_pw_data.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>

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
    void SetupPhon(const ModuleBase::Vector3<double> &k_d, const ModuleBase::Vector3<double> &q_d)
    {
        dfpt_test::setup_bases(pw_rho_, pw_wfc_, qlist_, data_, ucell_, k_d, q_d, 2, q_cart_);
        pert_.init(ucell_, &pw_rho_, &pw_wfc_, sf_);
        phon_.init(ucell_, &pw_rho_, &pert_);
    }
};

TEST_F(DFPTPhonSerialTest, IonIonAcousticSumRuleGamma)
{
    // a two-atom cell with different charges/masses breaks every symmetry:
    // the Gamma acoustic sum rule is then a razor for the Ewald balance
    // (G part + real part + Gaussian self term)
    dfpt_test::make_two_atom_cell(ucell_);
    ModuleBase::ComplexMatrix dyn(6, 6, true);
    phon_.ion_ion(ModuleBase::Vector3<double>(0.0, 0.0, 0.0), dyn);

    double max_elem = 0.0;
    for (int i = 0; i < 6; ++i)
    {
        for (int j = 0; j < 6; ++j)
        {
            max_elem = std::max(max_elem, std::abs(dyn(i, j)));
        }
    }
    ASSERT_GT(max_elem, 1.0e-6);
    // acoustic sum rule for the mass-scaled matrix D = Phi/sqrt(M_i M_j):
    // sum_j Phi(i,j) = 0  =>  sum_j sqrt(M_j) D(i,j) = 0 for every row i
    double sqrtm[2] = {std::sqrt(12.0), std::sqrt(4.0)};
    for (int i = 0; i < 6; ++i)
    {
        std::complex<double> rowsum(0.0, 0.0);
        for (int j = 0; j < 6; ++j)
        {
            rowsum += sqrtm[j / 3] * dyn(i, j);
        }
        EXPECT_LT(std::abs(rowsum), 1.0e-6 * max_elem) << "row " << i << " sum " << std::abs(rowsum);
    }
    // Hermitian
    for (int i = 0; i < 6; ++i)
    {
        for (int j = i + 1; j < 6; ++j)
        {
            EXPECT_NEAR(std::abs(dyn(i, j) - std::conj(dyn(j, i))), 0.0, 1.0e-10 * max_elem);
        }
    }
}

TEST_F(DFPTPhonSerialTest, IonIonGammaAcousticZeroModes)
{
    // same two-atom cell: three acoustic eigenvalues vanish at Gamma
    dfpt_test::make_two_atom_cell(ucell_);
    ModuleBase::ComplexMatrix dyn(6, 6, true);
    phon_.ion_ion(ModuleBase::Vector3<double>(0.0, 0.0, 0.0), dyn);
    data_.set_dynmat(0, dyn);
    phon_.diagonalize(0, data_);
    const std::vector<double> freq = data_.get_phon_freq(0);
    ASSERT_EQ(freq.size(), 6u);
    // three acoustic modes vanish; the frequencies come back in signed
    // ascending order, and a net-charged cell can push optical modes
    // negative (they then sort before the acoustic triple), so identify the
    // acoustic modes by magnitude
    std::vector<double> mag(freq);
    for (int i = 0; i < 6; ++i)
    {
        mag[i] = std::abs(freq[i]);
    }
    std::sort(mag.begin(), mag.end());
    for (int i = 0; i < 3; ++i)
    {
        EXPECT_LT(mag[i], 5.0) << "acoustic mode " << i; // cm^-1
    }
}

TEST_F(DFPTPhonSerialTest, IonIonGenericQVsDirectSum)
{
    // two-atom cell at an incommensurate q: the Ewald result must agree
    // with a direct (unscreened) dipole-Hessian lattice sum, whose shell
    // oscillation e^{i q l} makes it convergent
    dfpt_test::make_two_atom_cell(ucell_);
    const ModuleBase::Vector3<double> tau1(0.0, 0.0, 0.0);
    const ModuleBase::Vector3<double> tau2(0.25, 0.31, 0.17);
    const double z[2] = {4.0, 2.0};
    const double m[2] = {12.0, 4.0};
    const int nshell = 8; // lattice-vector cutoff in cells

    ModuleBase::ComplexMatrix dyn(6, 6, true);
    phon_.ion_ion(q_d_, dyn);

    // direct reference (structure validated against standalone Ewald sums):
    // off-diagonal (ia != ib):
    //   D_ab = -ZaZb e2/sqrt(MaMb) sum_l h0(R) e^{i2pi q.l}
    // diagonal: the on-site cross pairs are phase-free (both derivatives act
    // on tau_a in cell 0) while the same-atom images carry the phase
    // difference:
    //   D_aa = sum_{b != a} ZaZb e2/Ma sum_l h0(R)
    //        + Za^2 e2/Ma sum_{l != 0} h0(L) (1 - e^{i2pi q.l})
    // (h0(L) is the bare 1/R Hessian of the pure lattice; its l = 0 term is
    // killed by 1 - e^{i2pi q.0} = 0). The production Ewald drops the
    // tau-independent G = 0 constant (4pi/Omega)/3 on the on-site diagonal,
    // an ASR-preserving convention difference of ~2.5e-3 here, well inside
    // the tolerance below.
    ModuleBase::ComplexMatrix ref(6, 6, true);
    for (int ia = 0; ia < 2; ++ia)
    {
        for (int ib = 0; ib < 2; ++ib)
        {
            const bool self = (ib == ia);
            const ModuleBase::Vector3<double> dt = (ib == 0 ? tau1 : tau2) - (ia == 0 ? tau1 : tau2);
            for (int n1 = -nshell; n1 <= nshell; ++n1)
            {
                for (int n2 = -nshell; n2 <= nshell; ++n2)
                {
                    for (int n3 = -nshell; n3 <= nshell; ++n3)
                    {
                        if (self && n1 == 0 && n2 == 0 && n3 == 0)
                        {
                            continue;
                        }
                        const ModuleBase::Vector3<double> r(
                            (n1 * dfpt_test::cubic_edge + dt.x) * dfpt_test::lat0,
                            (n2 * dfpt_test::cubic_edge + dt.y) * dfpt_test::lat0,
                            (n3 * dfpt_test::cubic_edge + dt.z) * dfpt_test::lat0);
                        const double r2 = r * r;
                        const double r5 = r2 * r2 * std::sqrt(r2);
                        const double ph = ModuleBase::TWO_PI * (q_d_.x * n1 + q_d_.y * n2 + q_d_.z * n3);
                        const std::complex<double> phase(std::cos(ph), std::sin(ph));
                        const double pref = -z[ia] * z[ib] * ModuleBase::e2 / std::sqrt(m[ia] * m[ib]);
                        for (int da = 0; da < 3; ++da)
                        {
                            for (int db = 0; db < 3; ++db)
                            {
                                const double delta = (da == db) ? 1.0 : 0.0;
                                const double h0 = (3.0 * r[da] * r[db] - delta * r2) / r5;
                                if (self)
                                {
                                    ref(3 * ia + da, 3 * ia + db) +=
                                        z[ia] * z[ia] * ModuleBase::e2 / m[ia] * h0 * (1.0 - phase);
                                }
                                else
                                {
                                    ref(3 * ia + da, 3 * ib + db) += pref * h0 * phase;
                                    ref(3 * ia + da, 3 * ia + db) -= pref * std::sqrt(m[ib] / m[ia]) * h0;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    double max_ref = 0.0;
    for (int i = 0; i < 6; ++i)
    {
        for (int j = 0; j < 6; ++j)
        {
            max_ref = std::max(max_ref, std::abs(ref(i, j)));
        }
    }
    ASSERT_GT(max_ref, 1.0e-3);
    for (int i = 0; i < 6; ++i)
    {
        for (int j = 0; j < 6; ++j)
        {
            EXPECT_LT(std::abs(dyn(i, j) - ref(i, j)), 2.0e-3 * max_ref)
                << "(" << i << "," << j << ") ewald " << dyn(i, j) << " ref " << ref(i, j);
        }
    }
}

TEST_F(DFPTPhonSerialTest, IonIonFccSiQHalfEwaldConsistency)
{
    // fcc Si cell (2 atoms, quarter-diagonal basis) at the TR-invariant
    // q=(0.5,0,0): internal-consistency checks of the Ewald pair machinery
    // -- Hermiticity, the Gamma anchor value, alpha independence (two
    // pinned alphas give the same matrix), and determinism (two identical
    // calls give bit-equal matrices).
    // (A bare direct-sum reference for the off-diagonal blocks at this q
    // was attempted but its unoscillated cross-neighbor diagonal terms
    // are conditionally convergent and convention-dependent; the
    // element-exact agreement with a verbatim QE d2ionq transcription
    // was verified externally instead.)
    const double a_edge = 10.0;                // fcc cubic edge in lat0 units
    const double half_edge = 0.5 * a_edge;     // face-diagonal component of the fcc primitive vectors
    const double quarter_edge = 0.25 * a_edge; // basis-atom offset along the cubic diagonal
    ucell_.latvec =
        ModuleBase::Matrix3(0.0, half_edge, half_edge, half_edge, 0.0, half_edge, half_edge, half_edge, 0.0);
    ucell_.GT = ucell_.latvec.Inverse();
    ucell_.G = ucell_.GT.Transpose();
    ucell_.omega = a_edge * a_edge * a_edge * dfpt_test::lat0 * dfpt_test::lat0 * dfpt_test::lat0 / 4.0;
    dfpt_test::make_two_atom_cell(ucell_);
    ucell_.atoms[0].tau[0] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
    ucell_.atoms[1].tau[0] = ModuleBase::Vector3<double>(quarter_edge, quarter_edge, quarter_edge);
    // same species (Si-like)
    ucell_.atoms[1].ncpp.zv = 4.0;
    ucell_.atoms[1].mass = 12.0;
    // rebuild the bases for the new lattice (the SetUp bases used the
    // simple-cubic default cell)
    SetupPhon(k_d_, ModuleBase::Vector3<double>(0.5, 0.0, 0.0));

    const ModuleBase::Vector3<double> q_f(0.5, 0.0, 0.0);
    ModuleBase::ComplexMatrix dyn(6, 6, true);
    phon_.ion_ion(q_f, dyn);
    double max_elem = 0.0;
    for (int i = 0; i < 6; ++i)
    {
        for (int j = 0; j < 6; ++j)
        {
            max_elem = std::max(max_elem, std::abs(dyn(i, j)));
        }
    }
    ASSERT_GT(max_elem, 1.0e-6);
    // Hermitian
    for (int i = 0; i < 6; ++i)
    {
        for (int j = i + 1; j < 6; ++j)
        {
            EXPECT_NEAR(std::abs(dyn(i, j) - std::conj(dyn(j, i))), 0.0, 1.0e-12 * max_elem);
        }
    }
    // (no acoustic-sum-rule check here: ASR holds at Gamma only; the
    // dedicated IonIonAcousticSumRuleGamma test covers it)
    // Gamma anchor: the (1,2) block at Gamma is isotropic diagonal for this
    // cell; pin its value (validated against a verbatim QE d2ionq oracle)
    ModuleBase::ComplexMatrix dynG(6, 6, true);
    phon_.ion_ion(ModuleBase::Vector3<double>(0.0, 0.0, 0.0), dynG);
    EXPECT_NEAR(dynG(0, 3).real(), -0.0066210, 5.0e-6);
    EXPECT_NEAR(dynG(0, 4).real(), 0.0, 1.0e-8);
    // at q: the (1,2) block diagonal elements vanish and the off-diagonal
    // pattern is uniform (values pinned against the QE oracle)
    EXPECT_NEAR(std::abs(dyn(0, 3)), 0.0, 1.0e-8 * max_elem);
    EXPECT_NEAR(dyn(0, 4).real(), -0.0110598, 5.0e-6);
    EXPECT_NEAR(dyn(1, 0).real(), -0.0057142, 5.0e-6);
}

TEST_F(DFPTPhonSerialTest, IonIonAlphaIndependenceAndDeterminism)
{
    // alpha independence: pinning two different screening alphas (each
    // with its matching real-space radius) must give the same force
    // constants -- the Ewald split is exact for any converged alpha.
    // determinism: two identical calls must give bit-equal matrices.
    // Same fcc Si cell as IonIonFccSiQHalfEwaldConsistency.
    const double a_edge = 10.0;
    const double half_edge = 0.5 * a_edge;
    const double quarter_edge = 0.25 * a_edge;
    ucell_.latvec =
        ModuleBase::Matrix3(0.0, half_edge, half_edge, half_edge, 0.0, half_edge, half_edge, half_edge, 0.0);
    ucell_.GT = ucell_.latvec.Inverse();
    ucell_.G = ucell_.GT.Transpose();
    ucell_.omega = a_edge * a_edge * a_edge * dfpt_test::lat0 * dfpt_test::lat0 * dfpt_test::lat0 / 4.0;
    dfpt_test::make_two_atom_cell(ucell_);
    ucell_.atoms[0].tau[0] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
    ucell_.atoms[1].tau[0] = ModuleBase::Vector3<double>(quarter_edge, quarter_edge, quarter_edge);
    ucell_.atoms[1].ncpp.zv = 4.0;
    ucell_.atoms[1].mass = 12.0;
    SetupPhon(k_d_, ModuleBase::Vector3<double>(0.5, 0.0, 0.0));

    const ModuleBase::Vector3<double> q_f(0.5, 0.0, 0.0);
    ModuleBase::ComplexMatrix dyn_ref(6, 6, true);
    phon_.ion_ion(q_f, dyn_ref);

    // alpha independence: pin two distinct alphas (the rcut follows the
    // alpha through the ion_ion internals). Both values sit at or below
    // the internal search solution, so the fixed test-grid G ball still
    // bounds the reciprocal tail at the 1e-9 level; a larger alpha would
    // need a correspondingly larger ecutwfc to keep that tail bounded.
    ModuleBase::ComplexMatrix dyn_a(6, 6, true);
    phon_.set_ewald_alpha_override(0.15);
    phon_.ion_ion(q_f, dyn_a);
    ModuleBase::ComplexMatrix dyn_b(6, 6, true);
    phon_.set_ewald_alpha_override(0.30);
    phon_.ion_ion(q_f, dyn_b);
    for (int i = 0; i < 6; ++i)
    {
        for (int j = 0; j < 6; ++j)
        {
            EXPECT_NEAR(std::abs(dyn_a(i, j) - dyn_ref(i, j)), 0.0, 2.0e-6 * std::max(1.0, std::abs(dyn_ref(i, j))));
            EXPECT_NEAR(std::abs(dyn_b(i, j) - dyn_ref(i, j)), 0.0, 2.0e-6 * std::max(1.0, std::abs(dyn_ref(i, j))));
        }
    }

    // determinism: repeat the reference call, bit-equal
    ModuleBase::ComplexMatrix dyn_rep(6, 6, true);
    phon_.ion_ion(q_f, dyn_rep);
    for (int i = 0; i < 6; ++i)
    {
        for (int j = 0; j < 6; ++j)
        {
            EXPECT_DOUBLE_EQ(dyn_rep(i, j).real(), dyn_ref(i, j).real());
            EXPECT_DOUBLE_EQ(dyn_rep(i, j).imag(), dyn_ref(i, j).imag());
        }
    }
}

TEST_F(DFPTPhonSerialTest, EwaldAlphaSearchLowCutoffContinuesBelowDecade)
{
    // the multiplicative shrinkage must continue below alpha = 0.1: for
    // charge 8 and ggecut*tpiba2 = 8 the 1e-9 tail bound is first met
    // around alpha = 0.08, which a fixed-decrement search never reaches
    const double alpha = ModuleDFPT::ewald_alpha_search(8.0, 8.0, 1.0);
    EXPECT_GT(alpha, 0.0);
    EXPECT_LT(alpha, 0.1);
}
