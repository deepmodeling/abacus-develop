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

#define private public
#include "source_cell/atom_pseudo.h"
#include "source_cell/atom_spec.h"
#include "source_cell/pseudo.h"
#include "source_cell/qlist.h"
#include "source_cell/unitcell.h"
#include "source_cell/magnetism.h"
#include "source_pw/module_pwdft/stru_fac.h"
#include "source_pw/module_dfpt/dfpt_pert.h"
#include "source_pw/module_dfpt/dfpt_phon.h"
#undef private

#include "source_base/complexmatrix.h"
#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "source_base/vector3.h"
#include "source_lcao/module_dftu/dftu.h"
#include "source_psi/psi.h"

// test-support ctor/dtor stubs (see dfpt_pert_serial_test.cpp)
pseudo::pseudo()
{
}
pseudo::~pseudo()
{
}
Atom::Atom()
{
}
Atom::~Atom()
{
}
Atom_pseudo::Atom_pseudo()
{
}
Atom_pseudo::~Atom_pseudo()
{
}
SepPot::SepPot()
{
}
SepPot::~SepPot()
{
}
Sep_Cell::Sep_Cell() noexcept
{
}
Sep_Cell::~Sep_Cell() noexcept
{
}
UnitCell::UnitCell()
{
}
UnitCell::~UnitCell()
{
}
Magnetism::Magnetism()
{
}
Magnetism::~Magnetism()
{
}
Structure_Factor::Structure_Factor()
{
}
Structure_Factor::~Structure_Factor()
{
}

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

class DFPTPhonSerialTest : public testing::Test
{
  protected:
    const double lat0_ = 1.8897261254578281;
    const double ecutwfc_ = 2.5;
    const double rho_mult_ = 9.0;
    // cubic cell in lat0 units
    const double a_ = 10.0;

    ModuleBase::Matrix3 latvec_;
    UnitCell ucell_;
    ModulePW::PW_Basis pw_rho_;
    ModulePW::PW_Basis_K pw_wfc_;
    Structure_Factor sf_;
    ModuleDFPT::DFPT_Pert pert_;
    ModuleDFPT::DFPT_Phon phon_;
    ModuleCell::QList qlist_;
    ModuleDFPT::DFPT_PW_Data data_;

    const ModuleBase::Vector3<double> q_d_{0.13, 0.0, 0.07};
    const ModuleBase::Vector3<double> k_d_{-0.13, 0.0, -0.07};
    ModuleBase::Vector3<double> q_cart_;
    const ModuleBase::Vector3<double> tau_{1.1, 2.3, 0.7};

    void SetUp() override
    {
        latvec_ = ModuleBase::Matrix3(a_, 0.0, 0.0, 0.0, a_, 0.0, 0.0, 0.0, a_);
        ucell_.ntype = 1;
        ucell_.nat = 1;
        ucell_.atoms = new Atom[1];
        ucell_.atoms[0].na = 1;
        ucell_.atoms[0].tau.resize(1);
        ucell_.atoms[0].tau[0] = tau_;
        ucell_.latvec = latvec_;
        ucell_.GT = latvec_.Inverse();
        ucell_.G = ucell_.GT.Transpose();
        ucell_.lat0 = lat0_;
        ucell_.tpiba = ModuleBase::TWO_PI / lat0_;
        ucell_.tpiba2 = ucell_.tpiba * ucell_.tpiba;
        ucell_.omega = a_ * a_ * a_ * lat0_ * lat0_ * lat0_;
        ucell_.iat2it = new int[1];
        ucell_.iat2ia = new int[1];
        ucell_.iat2it[0] = 0;
        ucell_.iat2ia[0] = 0;
        MakeCoulombAtom();

        SetupBases(k_d_, q_d_);
    }

    // (re)initialize the bases and module wiring for a given (k, q) pair;
    // SetUp uses the default (k_d_, q_d_) fixture values
    void SetupBases(const ModuleBase::Vector3<double>& k_d,
                    const ModuleBase::Vector3<double>& q_d)
    {
        pw_rho_.initgrids(lat0_, latvec_, rho_mult_ * ecutwfc_);
        pw_rho_.initparameters(false, rho_mult_ * ecutwfc_);
        pw_rho_.fft_bundle.initfftmode(0);
        pw_rho_.setuptransform();
        pw_rho_.collect_local_pw();

        const ModuleBase::Vector3<double> klist[1] = {k_d};
        pw_wfc_.initgrids(lat0_, latvec_, pw_rho_.nx, pw_rho_.ny, pw_rho_.nz);
        pw_wfc_.initparameters(false, ecutwfc_, 1, klist);
        pw_wfc_.fft_bundle.initfftmode(0);
        pw_wfc_.setuptransform();
        pw_wfc_.collect_local_pw();

        qlist_.nkstot = 1;
        qlist_.kvec_d.clear();
        qlist_.kvec_d.push_back(q_d);
        q_cart_ = q_d * ucell_.G;

        data_.init(&qlist_, 1, 2, pw_wfc_.npwk_max, pw_rho_.nrxx, 1, 1, nullptr);
        pert_.init(ucell_, &pw_rho_, &pw_wfc_, sf_);
        phon_.init(ucell_, &pw_rho_, &pert_);
    }

    void TearDown() override
    {
        delete[] ucell_.atoms;
        ucell_.atoms = nullptr;
        delete[] ucell_.iat2it;
        ucell_.iat2it = nullptr;
        delete[] ucell_.iat2ia;
        ucell_.iat2ia = nullptr;
    }

    void MakeCoulombAtom()
    {
        Atom& at = ucell_.atoms[0];
        at.label = "C";
        at.coulomb_potential = true;
        at.ncpp.zv = 4.0;
        at.ncpp.tvanp = false;
        at.ncpp.has_so = false;
        at.ncpp.nbeta = 0;
        at.ncpp.nh = 0;
        at.ncpp.msh = 0;
        at.ncpp.kkbeta = 0;
        at.mass = 12.0;
    }

    double VlocCoulomb(double g2_bohr) const
    {
        return -ucell_.atoms[0].ncpp.zv * ModuleBase::e2 * ModuleBase::FOUR_PI / ucell_.omega / g2_bohr;
    }

    // reconfigure the cell as a two-atom Z=4/Z=2 crystal breaking all symmetry
    void MakeTwoAtomCell()
    {
        ucell_.ntype = 2;
        ucell_.nat = 2;
        delete[] ucell_.atoms;
        ucell_.atoms = new Atom[2];
        ucell_.atoms[0].na = 1;
        ucell_.atoms[1].na = 1;
        ucell_.atoms[0].tau.resize(1);
        ucell_.atoms[1].tau.resize(1);
        ucell_.atoms[0].tau[0] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
        ucell_.atoms[1].tau[0] = ModuleBase::Vector3<double>(0.25, 0.31, 0.17);
        for (int it = 0; it < 2; ++it)
        {
            Atom& at = ucell_.atoms[it];
            at.label = (it == 0) ? "A" : "B";
            at.coulomb_potential = true;
            at.ncpp.zv = (it == 0) ? 4.0 : 2.0;
            at.ncpp.tvanp = false;
            at.ncpp.has_so = false;
            at.ncpp.nbeta = 0;
            at.ncpp.nh = 0;
            at.mass = (it == 0) ? 12.0 : 4.0;
        }
        delete[] ucell_.iat2it;
        delete[] ucell_.iat2ia;
        ucell_.iat2it = new int[2];
        ucell_.iat2ia = new int[2];
        ucell_.iat2it[0] = 0;
        ucell_.iat2ia[0] = 0;
        ucell_.iat2it[1] = 1;
        ucell_.iat2ia[1] = 0;
    }

    // independent Ry/bohr^2/amu -> cm^-1 conversion used by diagonalize
    double RyBohr2AmuToCm1() const
    {
        const double amu_kg = 1.66053906660e-27; // CODATA amu in kg
        return std::sqrt(ModuleBase::RYDBERG_SI / amu_kg)
               / (0.529177210903e-10 * 2.0 * ModuleBase::PI * 2.99792458e10);
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
        EXPECT_LT(std::abs(rowsum), 1.0e-6 * max_elem)
            << "row " << i << " sum " << std::abs(rowsum);
    }
    // Hermitian
    for (int i = 0; i < 6; ++i)
    {
        for (int j = i + 1; j < 6; ++j)
        {
            EXPECT_NEAR(std::abs(dyn(i, j) - std::conj(dyn(j, i))),
                        0.0,
                        1.0e-10 * max_elem);
        }
    }
}

TEST_F(DFPTPhonSerialTest, IonIonGammaAcousticZeroModes)
{
    // same two-atom cell: three acoustic eigenvalues vanish at Gamma
    MakeTwoAtomCell();
    data_.set_dynmat(0, ModuleBase::ComplexMatrix(6, 6, true));
    ModuleBase::ComplexMatrix& dyn = data_.dynmat_[0];
    phon_.ion_ion(ModuleBase::Vector3<double>(0.0, 0.0, 0.0), dyn);
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
            const ModuleBase::Vector3<double> dt =
                (ib == 0 ? tau1 : tau2) - (ia == 0 ? tau1 : tau2);
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
                            (n1 * a_ + dt.x) * lat0_,
                            (n2 * a_ + dt.y) * lat0_,
                            (n3 * a_ + dt.z) * lat0_);
                        const double r2 = r * r;
                        const double r5 = r2 * r2 * std::sqrt(r2);
                        const double ph = ModuleBase::TWO_PI
                                          * (q_d_.x * n1 + q_d_.y * n2 + q_d_.z * n3);
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
                                        += z[ia] * z[ia] * ModuleBase::e2 / m[ia] * h0
                                           * (1.0 - phase);
                                }
                                else
                                {
                                    ref(3 * ia + da, 3 * ib + db) += pref * h0 * phase;
                                    ref(3 * ia + da, 3 * ia + db)
                                        -= pref * std::sqrt(m[ib] / m[ia]) * h0;
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
    const int npwk = pw_wfc_.npwk[0];
    psi::Psi<std::complex<double>> psi(1, 2, npwk, npwk, true);
    // the buffer is allocated uninitialized; zero it so only the components
    // set below are nonzero regardless of heap history from earlier tests
    psi.zero_out();
    // locate the G=0 plane wave in the k ball: getgpluskcar returns the
    // cartesian k+G (in 2pi/lat0 units), so look for k+G = k_cart, i.e. G = 0
    const ModuleBase::Vector3<double> k_cart = k_d_ * ucell_.G;
    int ig_zero = -1;
    for (int ig = 0; ig < npwk; ++ig)
    {
        const ModuleBase::Vector3<double> gk = pw_wfc_.getgpluskcar(0, ig);
        if (std::abs(gk.x - k_cart.x) < 1e-10 && std::abs(gk.y - k_cart.y) < 1e-10
            && std::abs(gk.z - k_cart.z) < 1e-10)
        {
            ig_zero = ig;
            break;
        }
    }
    ASSERT_GE(ig_zero, 0);
    psi(0, 0, ig_zero) = std::complex<double>(1.0, 0.0);
    ModuleBase::matrix wg(1, 2, true);
    wg(0, 0) = 2.0;
    wg(0, 1) = 0.0;

    // inject a known dpsi for displacement (atom 0, dir=1) on the k+q basis
    const int npwk_kq = [&]()
    {
        ModuleDFPT::DFPT_KQ_Basis kq;
        kq.init(&pw_wfc_, q_cart_, 0);
        return kq.get_npwk();
    }();
    std::vector<std::complex<double>> dpsi_inj(npwk_kq, std::complex<double>(0.0, 0.0));
    dpsi_inj[0] = std::complex<double>(0.3, 0.1);
    if (npwk_kq > 1)
    {
        dpsi_inj[1] = std::complex<double>(-0.2, 0.05);
    }
    data_.set_dpsi(0, 0, 0, dpsi_inj);

    phon_.accumulate_electron(0, 0, 1, psi, wg, data_);

    // expected: row 1 (atom 0, dir 1). The RHS on the k+q basis vector igl
    // carries the momentum w = Delta + q (Delta = G'' since k + q = 0 makes
    // every k+q basis vector a pure reciprocal-lattice harmonic G''):
    // RHS^a(G'') = -i tpiba w_a Vloc(w^2) e^{-i 2pi w.tau} (psi is a single
    // G'=0 plane wave and the Coulomb potential has no nonlocal part); the
    // GS structure-factor phase convention is exp(-i 2pi g.tau).
    ModuleDFPT::DFPT_KQ_Basis kq;
    kq.init(&pw_wfc_, q_cart_, 0);
    for (int adir = 0; adir < 3; ++adir)
    {
        std::complex<double> expect_cross(0.0, 0.0);
        for (int igl = 0; igl < kq.get_npwk(); ++igl)
        {
            const ModuleBase::Vector3<double> g = kq.get_gpluskq(igl); // = G''
            const ModuleBase::Vector3<double> w = g + q_cart_;          // Delta + q
            const double w2 = w * w;
            if (w2 < 1.0e-12)
            {
                continue; // Delta + q = 0 component dropped by dVloc
            }
            const double arg = -ModuleBase::TWO_PI * (w * tau_);
            const std::complex<double> rhs = std::complex<double>(0.0, -1.0)
                                             * (ucell_.tpiba * w[adir])
                                             * VlocCoulomb(w2 * ucell_.tpiba2)
                                             * std::complex<double>(std::cos(arg), std::sin(arg));
            expect_cross += std::conj(dpsi_inj[igl]) * rhs;
        }
        // the same-atom anharmonic d2 term: at this q = (0.13, 0, 0.07) the
        // second-order potential carries 2q = (0.26, 0, 0.14), which is NOT
        // a reciprocal vector: the same-k expectation is momentum-forbidden
        // and the production gate skips the whole term (zero contribution;
        // see the commensurate test below for the gate-on branch)
        // Hermitian 2n+1 accumulation: the row element receives
        // wg*<dpsi|dV^a|psi> once per off-diagonal column; the diagonal
        // column additionally gets its own conjugate (2 Re)
        std::complex<double> expect = wg(0, 0) * expect_cross;
        if (adir == 1)
        {
            expect = 2.0 * expect.real();
        }
        expect /= ucell_.atoms[0].mass;
        EXPECT_NEAR(std::abs(phon_.dynmat_accum_(1, adir) - expect),
                    0.0,
                    1.0e-7 * (1.0 + std::abs(expect)))
            << "adir " << adir << " got " << phon_.dynmat_accum_(1, adir)
            << " expect " << expect;
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
    const int npwk = pw_wfc_.npwk[0];
    psi::Psi<std::complex<double>> psi(1, 2, npwk, npwk, true);
    psi.zero_out();
    const ModuleBase::Vector3<double> k_cart = k_d_ * ucell_.G;
    int ig_zero = -1;
    for (int ig = 0; ig < npwk; ++ig)
    {
        const ModuleBase::Vector3<double> gk = pw_wfc_.getgpluskcar(0, ig);
        if (std::abs(gk.x - k_cart.x) < 1e-10 && std::abs(gk.y - k_cart.y) < 1e-10
            && std::abs(gk.z - k_cart.z) < 1e-10)
        {
            ig_zero = ig;
            break;
        }
    }
    ASSERT_GE(ig_zero, 0);
    psi(0, 0, ig_zero) = std::complex<double>(1.0, 0.0);
    ModuleBase::matrix wg(1, 2, true);
    wg(0, 0) = 2.0;
    wg(0, 1) = 0.0;

    ModuleDFPT::DFPT_KQ_Basis kq;
    kq.init(&pw_wfc_, q_cart_, 0);
    std::vector<std::complex<double>> dpsi_inj(kq.get_npwk(),
                                               std::complex<double>(0.0, 0.0));
    dpsi_inj[0] = std::complex<double>(0.25, -0.15);
    if (kq.get_npwk() > 2)
    {
        dpsi_inj[2] = std::complex<double>(0.4, 0.2);
    }
    data_.set_dpsi(0, 0, 0, dpsi_inj);

    phon_.accumulate_electron(0, 0, 0, psi, wg, data_);

    for (int adir = 0; adir < 3; ++adir)
    {
        std::complex<double> expect_cross(0.0, 0.0);
        for (int igl = 0; igl < kq.get_npwk(); ++igl)
        {
            const ModuleBase::Vector3<double> g = kq.get_gpluskq(igl);
            const ModuleBase::Vector3<double> w = g + q_cart_;
            const double w2 = w * w;
            if (w2 < 1.0e-12)
            {
                continue;
            }
            const double arg = -ModuleBase::TWO_PI * (w * tau_);
            const std::complex<double> rhs = std::complex<double>(0.0, -1.0)
                                             * (ucell_.tpiba * w[adir])
                                             * VlocCoulomb(w2 * ucell_.tpiba2)
                                             * std::complex<double>(std::cos(arg),
                                                                    std::sin(arg));
            expect_cross += std::conj(dpsi_inj[igl]) * rhs;
        }
        std::complex<double> expect = wg(0, 0) * expect_cross;
        if (adir == 0)
        {
            expect = 2.0 * expect.real();
        }
        expect /= ucell_.atoms[0].mass;
        EXPECT_NEAR(std::abs(phon_.dynmat_accum_(0, adir) - expect),
                    0.0,
                    1.0e-7 * (1.0 + std::abs(expect)))
            << "adir " << adir << " got " << phon_.dynmat_accum_(0, adir)
            << " expect " << expect;
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
    SetupBases(k_d, q_d);

    const int npwk = pw_wfc_.npwk[0];
    psi::Psi<std::complex<double>> psi(1, 2, npwk, npwk, true);
    psi.zero_out();
    const ModuleBase::Vector3<double> k_cart = k_d * ucell_.G;
    // three plane-wave components of psi at G' = 0, (0,1,0), (0,0,1) (all
    // inside the ecutwfc ball at this k): the d2 expectation lives on the
    // pairwise differences of |psi|^2 (the G=0 diagonal difference hits the
    // w=0 skip of the kernel); the (0,-1,1) difference makes the mixed
    // component K_{2,1} nonzero as well
    const int ncomp = 3;
    const ModuleBase::Vector3<double> gfrac[ncomp]
        = {ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
           ModuleBase::Vector3<double>(0.0, 1.0, 0.0),
           ModuleBase::Vector3<double>(0.0, 0.0, 1.0)};
    const std::complex<double> ccoef[ncomp] = {std::complex<double>(1.0, 0.0),
                                               std::complex<double>(0.6, -0.3),
                                               std::complex<double>(-0.4, 0.25)};
    ModuleBase::Vector3<double> gcart[ncomp];
    int ig_of[ncomp] = {-1, -1, -1};
    for (int ic = 0; ic < ncomp; ++ic)
    {
        gcart[ic] = gfrac[ic] * ucell_.G;
    }
    for (int ig = 0; ig < npwk; ++ig)
    {
        const ModuleBase::Vector3<double> gprim
            = pw_wfc_.getgpluskcar(0, ig) - k_cart;
        for (int ic = 0; ic < ncomp; ++ic)
        {
            if (std::abs(gprim.x - gcart[ic].x) < 1e-10
                && std::abs(gprim.y - gcart[ic].y) < 1e-10
                && std::abs(gprim.z - gcart[ic].z) < 1e-10)
            {
                ig_of[ic] = ig;
            }
        }
    }
    for (int ic = 0; ic < ncomp; ++ic)
    {
        ASSERT_GE(ig_of[ic], 0);
        psi(0, 0, ig_of[ic]) = ccoef[ic];
    }
    ModuleBase::matrix wg(1, 2, true);
    wg(0, 0) = 2.0;
    wg(0, 1) = 0.0;

    // injected dpsi on the k+q = 0 ball (arbitrary coefficients)
    ModuleDFPT::DFPT_KQ_Basis kq;
    kq.init(&pw_wfc_, q_cart_, 0);
    const int npwk_kq = kq.get_npwk();
    std::vector<std::complex<double>> dpsi_inj(npwk_kq, std::complex<double>(0.0, 0.0));
    dpsi_inj[0] = std::complex<double>(0.3, 0.1);
    if (npwk_kq > 1)
    {
        dpsi_inj[1] = std::complex<double>(-0.2, 0.05);
    }
    data_.set_dpsi(0, 0, 0, dpsi_inj);

    phon_.accumulate_electron(0, 0, 1, psi, wg, data_);

    for (int adir = 0; adir < 3; ++adir)
    {
        // cross term: RHS(G'') = sum_i c_i (-i) tpiba w_{i,a} Vloc(|w_i|^2)
        // e^{-i2pi w_i.tau} with w_i = G'' - G'_i + q
        std::complex<double> expect_cross(0.0, 0.0);
        for (int igl = 0; igl < npwk_kq; ++igl)
        {
            const ModuleBase::Vector3<double> gpp = kq.get_gpluskq(igl); // = G''
            for (int ic = 0; ic < ncomp; ++ic)
            {
                const ModuleBase::Vector3<double> w = gpp - gcart[ic] + q_cart_;
                const double w2 = w * w;
                if (w2 < 1.0e-12)
                {
                    continue;
                }
                const double arg = -ModuleBase::TWO_PI * (w * tau_);
                const std::complex<double> rhs = std::complex<double>(0.0, -1.0)
                                                 * (ucell_.tpiba * w[adir])
                                                 * VlocCoulomb(w2 * ucell_.tpiba2)
                                                 * std::complex<double>(std::cos(arg),
                                                                        std::sin(arg));
                expect_cross += ccoef[ic] * std::conj(dpsi_inj[igl]) * rhs;
            }
        }
        std::complex<double> expect = wg(0, 0) * expect_cross;
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
            for (int i = 0; i < ncomp; ++i)
            {
                for (int j = 0; j < ncomp; ++j)
                {
                    const ModuleBase::Vector3<double> g = gcart[i] - gcart[j];
                    const double g2 = g * g;
                    if (g2 < 1.0e-12)
                    {
                        continue;
                    }
                    const double arg = -ModuleBase::TWO_PI * (g * tau_);
                    const std::complex<double> kterm
                        = -(ucell_.tpiba * g[adir]) * (ucell_.tpiba * g[1])
                          * VlocCoulomb(g2 * ucell_.tpiba2)
                          * std::complex<double>(std::cos(arg), std::sin(arg));
                    d2elem += std::conj(ccoef[i]) * ccoef[j] * kterm;
                }
            }
            expect += wg(0, 0) * d2elem;
        }
        expect /= ucell_.atoms[0].mass;
        EXPECT_NEAR(std::abs(phon_.dynmat_accum_(1, adir) - expect),
                    0.0,
                    1.0e-7 * (1.0 + std::abs(expect)))
            << "adir " << adir << " got " << phon_.dynmat_accum_(1, adir)
            << " expect " << expect;
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
    auto block = [&expect](double a, double b, std::complex<double> c)
    {
        const double mid = 0.5 * (a + b);
        const double rad = std::sqrt(std::pow(0.5 * (a - b), 2) + std::norm(c));
        expect.push_back(mid + rad);
        expect.push_back(mid - rad);
    };
    block(lam[0], lam[1], dyn(0, 1)); // coupled pair
    block(lam[2], lam[3], dyn(2, 3)); // coupled pair
    expect.push_back(lam[4]); // untouched diagonal
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

// ---------------------------------------------------------------------------
// add_loto / check_sum_rule
// ---------------------------------------------------------------------------

TEST_F(DFPTPhonSerialTest, AddLotoIsotropicClosedForm)
{
    // isotropic eps_inf = 3, Born charges Z*_1 = 1, Z*_2 = 2, masses 12/4
    ModuleBase::ComplexMatrix dyn0(6, 6, true);
    data_.set_dynmat(0, dyn0);
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
    // temporarily make the cell two-atom for mass lookup consistency
    ucell_.ntype = 2;
    ucell_.nat = 2;
    delete[] ucell_.atoms;
    ucell_.atoms = new Atom[2];
    ucell_.atoms[0].na = 1;
    ucell_.atoms[1].na = 1;
    ucell_.atoms[0].mass = 12.0;
    ucell_.atoms[1].mass = 4.0;
    delete[] ucell_.iat2it;
    delete[] ucell_.iat2ia;
    ucell_.iat2it = new int[2];
    ucell_.iat2ia = new int[2];
    ucell_.iat2it[0] = 0;
    ucell_.iat2ia[0] = 0;
    ucell_.iat2it[1] = 1;
    ucell_.iat2ia[1] = 0;

    const ModuleBase::Vector3<double> qhat(1.0, 0.0, 0.0);
    phon_.add_loto(qhat, data_);

    // closed form: D_NAC(0x,1x) = 4pi e2/Omega * 1*2/(3) / sqrt(12*4)
    const double expect = ModuleBase::FOUR_PI * ModuleBase::e2 / ucell_.omega / 3.0
                          * 2.0 / std::sqrt(48.0);
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
