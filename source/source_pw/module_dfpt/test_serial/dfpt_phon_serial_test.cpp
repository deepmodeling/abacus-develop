#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>

// serial unit test of the DFPT dynamical matrix (C5): the Ewald ion-ion
// force constants, the electronic 2n+1 accumulation, the Hermitian
// eigensolver and the LO-TO term. Runs without __MPI on the shared FFT grid
// like the other DFPT serial tests.

#include "dfpt_serial_fixture.h"
#include "source_base/complexmatrix.h"
#include "source_pw/module_dfpt/dfpt_rho.h"
#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "source_base/vector3.h"
#include "source_cell/atom_pseudo.h"
#include "source_cell/atom_spec.h"
#include "source_cell/qlist.h"
#include "source_cell/unitcell.h"
#include "source_psi/psi.h"
#include "source_pw/module_dfpt/dfpt_kq_basis.h"
#include "source_pw/module_dfpt/dfpt_pert.h"
#include "source_pw/module_dfpt/dfpt_phon.h"
#include "source_cell/stru_fac.h"

// ctor/dtor stubs for the cell/spepot/stru_fac link closures live in the
// shared test/dfpt_test_mocks.cpp compiled into every DFPT test binary.

/************************************************
 *  serial unit test of DFPT_Phon (C5)
 ***********************************************/

/**
 * - Tested Functions:
 *   - ion_ion: the Ewald force constants satisfy the acoustic sum rule at
 *     Gamma to grid accuracy (this pins the sign and magnitude of the
 *     Gaussian self term), give three zero acoustic modes at Gamma, and are
 *     cross-checked against a direct (unscreened) lattice Hessian sum at a
 *     generic incommensurate q where the dipole sum is oscillation-screened.
 *   - accumulate_electron: with an injected dpsi, the cross term
 *     2 sum wg Re <dpsi | dV^a_ext | psi> is validated against an analytic
 *     convolution (psi a single plane wave, Coulomb local potential), and
 *     the same-atom anharmonic term <psi|d2V|psi> against the closed-form
 *     coefficient at Delta = 0 (|u|^2 = 1 keeps only the G=0 harmonic);
 *     the dpsi slot is restored after the accumulation.
 *   - diagonalize: signed frequencies of a known complex Hermitian matrix
 *     against an independently computed Ry/bohr^2/amu -> cm^-1 factor.
 *   - add_loto: isotropic dielectric/Born-charge LO-TO term against the
 *     closed-form matrix element.
 *   - check_sum_rule at Gamma.
 */

class DFPTPhonSerialTest : public DFPTSerialBase
{
  protected:
    Structure_Factor sf_;
    ModuleDFPT::DFPT_Pert pert_;
    ModuleDFPT::DFPT_Phon phon_;

    void SetUp() override
    {
        DFPTSerialBase::SetUp();
        SetupPhon(k_d_, q_d_);
    }

    // (re)initialize the bases and the pert/phon wiring for a given (k, q)
    // pair; SetUp uses the default fixture values
    void SetupPhon(const ModuleBase::Vector3<double>& k_d, const ModuleBase::Vector3<double>& q_d)
    {
        SetupBases(k_d, q_d, 2);
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
        MakeTwoAtomCell();
    }

    // band 0 occupied with wg = 2, band 1 unoccupied
    static ModuleBase::matrix MakeOccWeights()
    {
        ModuleBase::matrix wg(1, 2, true);
        wg(0, 0) = 2.0;
        wg(0, 1) = 0.0;
        return wg;
    }

    // psi: band 0 = single plane wave at G' = 0 (c = 1, occupied),
    // band 1 unoccupied. The buffer is allocated uninitialized; zero it so
    // only the component set below is nonzero regardless of heap history
    // from earlier tests. getgpluskcar returns the cartesian k+G, so the
    // G = 0 entry is the one with k+G = k_cart.
    psi::Psi<std::complex<double>> MakeSinglePlaneWavePsi(const ModuleBase::Vector3<double>& k_d)
    {
        const int npwk = pw_wfc_.npwk[0];
        psi::Psi<std::complex<double>> psi(1, 2, npwk, npwk, true);
        psi.zero_out();
        const ModuleBase::Vector3<double> k_cart = k_d * ucell_.G;
        for (int ig = 0; ig < npwk; ++ig)
        {
            const ModuleBase::Vector3<double> gk = pw_wfc_.getgpluskcar(0, ig);
            if (std::abs(gk.x - k_cart.x) < 1e-10 && std::abs(gk.y - k_cart.y) < 1e-10
                && std::abs(gk.z - k_cart.z) < 1e-10)
            {
                psi(0, 0, ig) = std::complex<double>(1.0, 0.0);
                break;
            }
        }
        return psi;
    }

    // analytic cross term sum_G'' conj(dpsi_G'') sum_i c_i
    // RHS^a(G''-G'_i): the first-order local Coulomb potential on the k+q
    // basis, RHS^a(w) = -i tpiba w_a Vloc(|w|^2) e^{-i 2pi w.tau} with
    // w = G'' - G' + q (GS structure-factor phase convention); the
    // Delta + q = 0 component is dropped by dVloc
    std::complex<double> AnalyticCrossTerm(const ModuleDFPT::DFPT_KQ_Basis& kq,
                                           const std::vector<std::complex<double>>& psi_coef,
                                           const std::vector<ModuleBase::Vector3<double>>& psi_gcart,
                                           const std::vector<std::complex<double>>& dpsi_inj,
                                           int adir) const
    {
        std::complex<double> cross(0.0, 0.0);
        for (int igl = 0; igl < kq.get_npwk(); ++igl)
        {
            const ModuleBase::Vector3<double> gpp = kq.get_gpluskq(igl);
            for (size_t ic = 0; ic < psi_coef.size(); ++ic)
            {
                // AnalyticDVloc returns 0 at w = 0 (dVloc drop)
                cross += psi_coef[ic] * std::conj(dpsi_inj[igl]) * AnalyticDVloc(adir, gpp - psi_gcart[ic] + q_cart_);
            }
        }
        return cross;
    }
};

// ---------------------------------------------------------------------------
// ion_ion
// ---------------------------------------------------------------------------

TEST_F(DFPTPhonSerialTest, IonIonAcousticSumRuleGamma)
{
    // a two-atom cell with different charges/masses breaks every symmetry:
    // the Gamma acoustic sum rule is then a razor for the Ewald balance
    // (G part + real part + Gaussian self term)
    MakeTwoAtomCell();
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
    MakeTwoAtomCell();
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
    MakeTwoAtomCell();
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
                        const ModuleBase::Vector3<double> r((n1 * a_ + dt.x) * lat0_,
                                                            (n2 * a_ + dt.y) * lat0_,
                                                            (n3 * a_ + dt.z) * lat0_);
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
                                    ref(3 * ia + da, 3 * ia + db)
                                        += z[ia] * z[ia] * ModuleBase::e2 / m[ia] * h0 * (1.0 - phase);
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

// ---------------------------------------------------------------------------
// accumulate_electron
// ---------------------------------------------------------------------------

TEST_F(DFPTPhonSerialTest, AccumulateElectronAnalyticContraction)
{
    // psi: band 0 = single plane wave at G'=0 (c=1, occupied, wg=2),
    // band 1 unoccupied. k = -q so the k+q basis vectors are plain G''.
    psi::Psi<std::complex<double>> psi = MakeSinglePlaneWavePsi(k_d_);
    ModuleBase::matrix wg = MakeOccWeights();

    // inject a known dpsi for displacement (atom 0, dir=1) on the k+q basis
    ModuleDFPT::DFPT_KQ_Basis kq;
    kq.init(&pw_wfc_, &pw_rho_, q_cart_, 0);
    std::vector<std::complex<double>> dpsi_inj(kq.get_npwk(), std::complex<double>(0.0, 0.0));
    dpsi_inj[0] = std::complex<double>(0.3, 0.1);
    if (kq.get_npwk() > 1)
    {
        dpsi_inj[1] = std::complex<double>(-0.2, 0.05);
    }
    data_.set_dpsi(0, 0, 0, dpsi_inj);

    phon_.accumulate_electron(0, 0, 1, psi, wg, data_, std::vector<std::complex<double>>());

    // expected: row 1 (atom 0, dir 1) of the Hermitian 2n+1 accumulation:
    // the row element receives wg*<dpsi|dV^a|psi> once per off-diagonal
    // column; the diagonal column additionally gets its own conjugate
    // (2 Re). The same-atom anharmonic d2 term is gate-skipped here: at
    // this q = (0.13, 0, 0.07) the second-order potential carries
    // 2q = (0.26, 0, 0.14), which is NOT a reciprocal vector, so the
    // same-k expectation is momentum-forbidden (see the commensurate
    // test below for the gate-on branch)
    const std::vector<ModuleBase::Vector3<double>> g0(1, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    for (int adir = 0; adir < 3; ++adir)
    {
        std::complex<double> expect
            = wg(0, 0) * AnalyticCrossTerm(kq, {std::complex<double>(1.0, 0.0)}, g0, dpsi_inj, adir);
        if (adir == 1)
        {
            expect = 2.0 * expect.real();
        }
        expect /= ucell_.atoms[0].mass;
        EXPECT_NEAR(std::abs(phon_.dynmat_accum()(1, adir) - expect), 0.0, 1.0e-7 * (1.0 + std::abs(expect)))
            << "adir " << adir << " got " << phon_.dynmat_accum()(1, adir) << " expect " << expect;
    }

    // the dpsi slot must be restored to the injected solution
    const std::vector<std::complex<double>> restored = data_.get_dpsi(0, 0, 0);
    ASSERT_EQ(restored.size(), dpsi_inj.size());
    for (size_t i = 0; i < dpsi_inj.size(); ++i)
    {
        EXPECT_DOUBLE_EQ(restored[i].real(), dpsi_inj[i].real());
        EXPECT_DOUBLE_EQ(restored[i].imag(), dpsi_inj[i].imag());
    }
}

TEST_F(DFPTPhonSerialTest, AccumulateElectronD2GateOffGenericQ)
{
    // row 0 of the same generic-q fixture: the same-atom d2 kernel would be
    // nonzero here under the ungated convention (the (0,0) element involves
    // q_x^2 != 0), so this row is a sharp probe that the 2q-reciprocal gate
    // really suppresses the momentum-forbidden term at a generic q
    psi::Psi<std::complex<double>> psi = MakeSinglePlaneWavePsi(k_d_);
    ModuleBase::matrix wg = MakeOccWeights();

    ModuleDFPT::DFPT_KQ_Basis kq;
    kq.init(&pw_wfc_, &pw_rho_, q_cart_, 0);
    std::vector<std::complex<double>> dpsi_inj(kq.get_npwk(), std::complex<double>(0.0, 0.0));
    dpsi_inj[0] = std::complex<double>(0.25, -0.15);
    if (kq.get_npwk() > 2)
    {
        dpsi_inj[2] = std::complex<double>(0.4, 0.2);
    }
    data_.set_dpsi(0, 0, 0, dpsi_inj);

    phon_.accumulate_electron(0, 0, 0, psi, wg, data_, std::vector<std::complex<double>>());

    const std::vector<ModuleBase::Vector3<double>> g0(1, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    for (int adir = 0; adir < 3; ++adir)
    {
        std::complex<double> expect
            = wg(0, 0) * AnalyticCrossTerm(kq, {std::complex<double>(1.0, 0.0)}, g0, dpsi_inj, adir);
        if (adir == 0)
        {
            expect = 2.0 * expect.real();
        }
        expect /= ucell_.atoms[0].mass;
        EXPECT_NEAR(std::abs(phon_.dynmat_accum()(0, adir) - expect), 0.0, 1.0e-7 * (1.0 + std::abs(expect)))
            << "adir " << adir << " got " << phon_.dynmat_accum()(0, adir) << " expect " << expect;
    }
}

TEST_F(DFPTPhonSerialTest, AccumulateElectronD2CommensurateQ)
{
    // k = (-0.5, 0, 0) and q = (0.5, 0, 0): the k+q ball is centered at 0
    // (pure G'' harmonics for the cross term) and 2q = (1, 0, 0) IS a
    // reciprocal vector, so the same-atom d2 gate passes: the second-order
    // local operator is lattice-periodic on the integer G set (q_eff = 0)
    // with kernel K_{da,db}(G) = -tpiba^2 G_da G_db Vloc(|G|^2) e^{-i2pi G.tau}
    const ModuleBase::Vector3<double> k_d(-0.5, 0.0, 0.0);
    const ModuleBase::Vector3<double> q_d(0.5, 0.0, 0.0);
    SetupPhon(k_d, q_d);

    const int npwk = pw_wfc_.npwk[0];
    psi::Psi<std::complex<double>> psi(1, 2, npwk, npwk, true);
    psi.zero_out();
    const ModuleBase::Vector3<double> k_cart = k_d * ucell_.G;
    // three plane-wave components of psi at G' = 0, (0,1,0), (0,0,1) (all
    // inside the ecutwfc ball at this k): the d2 expectation lives on the
    // pairwise differences of |psi|^2 (the G=0 diagonal difference hits the
    // w=0 skip of the kernel); the (0,-1,1) difference makes the mixed
    // component K_{2,1} nonzero as well
    const std::vector<ModuleBase::Vector3<double>> gfrac = {ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                                                            ModuleBase::Vector3<double>(0.0, 1.0, 0.0),
                                                            ModuleBase::Vector3<double>(0.0, 0.0, 1.0)};
    const std::vector<std::complex<double>> ccoef
        = {std::complex<double>(1.0, 0.0), std::complex<double>(0.6, -0.3), std::complex<double>(-0.4, 0.25)};
    const size_t ncomp = gfrac.size();
    std::vector<ModuleBase::Vector3<double>> gcart(ncomp);
    std::vector<int> ig_of(ncomp, -1);
    for (size_t ic = 0; ic < ncomp; ++ic)
    {
        gcart[ic] = gfrac[ic] * ucell_.G;
    }
    for (int ig = 0; ig < npwk; ++ig)
    {
        const ModuleBase::Vector3<double> gprim = pw_wfc_.getgpluskcar(0, ig) - k_cart;
        for (size_t ic = 0; ic < ncomp; ++ic)
        {
            if (std::abs(gprim.x - gcart[ic].x) < 1e-10 && std::abs(gprim.y - gcart[ic].y) < 1e-10
                && std::abs(gprim.z - gcart[ic].z) < 1e-10)
            {
                ig_of[ic] = ig;
            }
        }
    }
    for (size_t ic = 0; ic < ncomp; ++ic)
    {
        ASSERT_GE(ig_of[ic], 0);
        psi(0, 0, ig_of[ic]) = ccoef[ic];
    }
    const ModuleBase::matrix wg = MakeOccWeights();

    // injected dpsi on the k+q = 0 ball (arbitrary coefficients)
    ModuleDFPT::DFPT_KQ_Basis kq;
    kq.init(&pw_wfc_, &pw_rho_, q_cart_, 0);
    std::vector<std::complex<double>> dpsi_inj(kq.get_npwk(), std::complex<double>(0.0, 0.0));
    dpsi_inj[0] = std::complex<double>(0.3, 0.1);
    if (kq.get_npwk() > 1)
    {
        dpsi_inj[1] = std::complex<double>(-0.2, 0.05);
    }
    data_.set_dpsi(0, 0, 0, dpsi_inj);

    phon_.accumulate_electron(0, 0, 1, psi, wg, data_, std::vector<std::complex<double>>());

    for (int adir = 0; adir < 3; ++adir)
    {
        std::complex<double> expect = wg(0, 0) * AnalyticCrossTerm(kq, ccoef, gcart, dpsi_inj, adir);
        if (adir == 1)
        {
            expect = 2.0 * expect.real();
        }
        // d2 term: gate passes at this q; the columns cola >= rowb = 1
        // (adir 1 and 2) receive wg * <psi|K_{adir,1}|psi>. The |u|^2
        // harmonic at G'_j - G'_i probes the kernel at the NEGATIVE
        // harmonic, so the closed form runs over K(G'_i - G'_j) with
        // coefficient c_i* c_j (K(-g) = conj K(g), K(0) = 0); the diagonal
        // column adds it once (real), the off-diagonal once
        if (adir >= 1)
        {
            std::complex<double> d2elem(0.0, 0.0);
            for (size_t i = 0; i < ncomp; ++i)
            {
                for (size_t j = 0; j < ncomp; ++j)
                {
                    const ModuleBase::Vector3<double> g = gcart[i] - gcart[j];
                    const double g2 = g * g;
                    if (g2 < 1.0e-12)
                    {
                        continue;
                    }
                    const double arg = -ModuleBase::TWO_PI * (g * tau_);
                    const std::complex<double> kterm = -(ucell_.tpiba * g[adir]) * (ucell_.tpiba * g[1])
                                                       * VlocCoulomb(g2 * ucell_.tpiba2)
                                                       * std::complex<double>(std::cos(arg), std::sin(arg));
                    d2elem += std::conj(ccoef[i]) * ccoef[j] * kterm;
                }
            }
            expect += wg(0, 0) * d2elem;
        }
        expect /= ucell_.atoms[0].mass;
        EXPECT_NEAR(std::abs(phon_.dynmat_accum()(1, adir) - expect), 0.0, 1.0e-7 * (1.0 + std::abs(expect)))
            << "adir " << adir << " got " << phon_.dynmat_accum()(1, adir) << " expect " << expect;
    }
}

// ---------------------------------------------------------------------------
// diagonalize
// ---------------------------------------------------------------------------

TEST_F(DFPTPhonSerialTest, DiagonalizeKnownMatrix)
{
    // 2-atom layout so nat3 = 6. Hermitian blocks with known closed-form
    // spectra: [[a, c], [conj(c), b]] has eigenvalues (a+b)/2
    // +- sqrt(((a-b)/2)^2 + |c|^2)
    MakeTwoAtomCell();
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
    for (double& e: expect)
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

// ---------------------------------------------------------------------------
// add_loto / check_sum_rule
// ---------------------------------------------------------------------------

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

// ---------------------------------------------------------------------------
// B2: LO-TO direction via the data layer, corrected-frequency
// diagonalization, and the output format regression
// ---------------------------------------------------------------------------

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

// ---------------------------------------------------------------------------
// nlcc direct-term fill-semantics regression: with a synthetic provider set
// (all-atoms kernel table, ground-state XC vector, DFPT_Rho), the nlcc
// kernel accumulation must satisfy the QE-pinned invariants:
//   (1) the one-sided cross block (nlcc row, coreless column) receives
//       (1/2)(X_ab + conj(X_ba)) = (1/2) X_ab  [X_ba = 0: coreless row]
//   (2) every displacement slot's dvaux participates (the "dead row" bug
//       dropped the entire X(N, B) cross block when the table left
//       coreless atoms empty)
//   (3) the matrix stays Hermitian
// These pin the fill factor (the double-sided same-atom block would double
// under a both-order fill; a lower-triangle-only fill halves the one-sided
// cross blocks).
// ---------------------------------------------------------------------------
TEST_F(DFPTPhonSerialTest, NlccDirectTermsFillSemantics)
{
    MakeNlccTwoAtomCell();
    SetupPhon(k_d_, q_d_);

    ModuleDFPT::DFPT_Rho rho;
    rho.init({1, pw_rho_.nrxx, &pw_rho_, &pw_wfc_, G_, std::string("rho_plain"), 0.4, 1.0, 8});
    phon_.set_rho(&rho);

    const int nrxx = pw_rho_.nrxx;
    // synthetic kernel table: f_xc.drhoc stand-ins = distinct constant
    // complex fields per slot (values chosen to survive exact comparison)
    std::vector<std::vector<std::complex<double>>> table(6);
    for (int slot = 0; slot < 6; ++slot)
    {
        const double re = 0.1 + 0.05 * slot;
        const double im = 0.02 * slot - 0.05;
        table[slot].assign(nrxx, std::complex<double>(re, im));
    }
    phon_.set_dvxc_core_table(&table);

    // synthetic ground-state XC vector (only the diag term consumes it);
    // zero it so the diag term's only effect is absent (this test targets
    // the kernel fill)
    const std::vector<double> vxc_zero(nrxx, 0.0);
    phon_.set_vxc_r(&vxc_zero);

    psi::Psi<std::complex<double>> psi = MakeSinglePlaneWavePsi(k_d_);
    ModuleBase::matrix wg = MakeOccWeights();
    // NO dpsi injection: the bare cross term needs occupied dpsi entries
    // and skips empty slots, so the accumulated matrix isolates the nlcc
    // kernel contribution exactly

    phon_.accumulate_electron(0, 1, 0, psi, wg, data_, table[3]);

    const auto& m = phon_.dynmat_accum();
    // row (1, dir 0) = 3*1+0 = 3; column (atom 0, dir j) = j (coreless)
    // X_3j = sum_r table[3](r) * conj(drhoc(atom0, dir_j)(r)) = 0
    //   (atom 0 has no nlcc: drhoc = empty -> the kernel term's x_ab arm
    //    is zero for coreless COLUMNS)
    // x_ba arm = sum_r table[j](r) * conj(drhoc(atom1, dir0)(r)) != 0
    //   (row atom 1 has a core; the column's dvaux comes from the table)
    // So the (3, j) element = (1/2)(0 + conj(x_ba)) = (1/2) conj(x_ba).
    // The x_ba contraction runs over the actual Gaussian core field, so
    // verify structurally: the one-sided block is NONZERO (dead-row bug
    // gave exactly zero) and matches an independent recomputation.
    const double fac = ucell_.omega / nrxx;
    std::vector<std::complex<double>> drhoc_g;
    rho.compute_drhoc_g(ucell_, 1, 0, q_cart_, drhoc_g);
    ASSERT_EQ(static_cast<int>(drhoc_g.size()), pw_rho_.npw);
    std::vector<std::complex<double>> drhoc_r(nrxx);
    pw_rho_.recip2real(drhoc_g.data(), drhoc_r.data());
    const double m0 = std::sqrt(ucell_.atoms[0].mass);
    const double m1 = std::sqrt(ucell_.atoms[1].mass);
    for (int j = 0; j < 3; ++j)
    {
        std::complex<double> x_ba(0.0, 0.0);
        for (int ir = 0; ir < nrxx; ++ir)
        {
            x_ba += std::conj(table[j][ir]) * drhoc_r[ir];
        }
        x_ba *= fac;
        const std::complex<double> expect = 0.5 * std::conj(x_ba) / (m0 * m1);
        EXPECT_NEAR(std::abs(m(3, j) - expect), 0.0, 1.0e-10 * (1.0 + std::abs(expect)))
            << "j=" << j << " got " << m(3, j) << " expect " << expect;
    }

    // Hermiticity of the accumulated block so far: accumulate the column
    // displacement (atom 1, dir 0 already done as row; do (atom 1, dir 1)
    // to complete more of the same-atom block, then check the whole
    // matrix is Hermitian within the accumulated rows
    phon_.accumulate_electron(0, 1, 1, psi, wg, data_, table[4]);
    // the (3, 4) and (4, 3) elements must be conjugate pairs: from row 3's
    // pass, (3,4) = (1/2)(X_34 + conj(X_43)); from row 4's pass, (4,3) =
    // (1/2)(X_43 + conj(X_34)) = conj((3,4))
    const std::complex<double> e34 = m(3, 4);
    const std::complex<double> e43 = m(4, 3);
    EXPECT_NEAR(e43.real(), e34.real(), 1.0e-12);
    EXPECT_NEAR(e43.imag(), -e34.imag(), 1.0e-12);
}
