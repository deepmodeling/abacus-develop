// Serial unit tests of DFPT_Phon::accumulate_electron (dfpt_phon_elec.cpp):
// with an injected dpsi, the cross term 2 sum wg Re <dpsi | dV^a_ext | psi>
// is validated against an analytic convolution (psi a single plane wave,
// Coulomb local potential), and the same-atom anharmonic term
// <psi|d2V|psi> against the closed-form coefficient at Delta = 0 (|u|^2
// = 1 keeps only the G=0 harmonic); the dpsi slot is restored after the
// accumulation. Plus the nlcc kernel fill-semantics regression.

#include "dfpt_phon_serial_fixture.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <cmath>
#include <complex>
#include <vector>

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
