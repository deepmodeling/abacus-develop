// Serial unit tests of the DFPT_Rho response-mixing family (the mix_drho
// stack inside dfpt_rho.cpp; this file's name anticipates the Phase-2
// source split that moves the mixing stack into dfpt_rho_mix.cpp):
// plain first/second-step combinations and the residual formula, the
// Kerker preconditioned first step, the stiff-model stabilization, and
// the Modified-Broyden convergence and fresh-q reset regressions.

#include "dfpt_rho_serial_fixture.h"

#include "gtest/gtest.h"

using namespace dfpt_rho_test_rng;

TEST_F(DFPTRhoSerialTest, MixDrhoFirstStepIsScaledOutput)
{
    psi::Psi<std::complex<double>> psi;
    std::vector<std::vector<std::vector<std::complex<double>>>> dpsi;
    FillRandomStates(psi, dpsi);
    ModuleBase::matrix wg(1, nbands_);
    wg(0, 0) = 1.0;
    wg(0, 1) = 0.0;

    rho_.compute_drho(psi, wg, 0, data_);
    const std::vector<std::complex<double>> out = data_.get_drho_g(0, 0);
    rho_.mix_drho(0, data_);
    const std::vector<std::complex<double>> mixed = data_.get_drho_g(0, 0);
    ASSERT_EQ(mixed.size(), out.size());
    for (int ig = 0; ig < pw_rho_.npw; ++ig)
    {
        EXPECT_NEAR(mixed[ig].real(), 0.4 * out[ig].real(), 1.0e-12);
        EXPECT_NEAR(mixed[ig].imag(), 0.4 * out[ig].imag(), 1.0e-12);
    }
    // zero input: ||out - 0|| / ||out|| == 1 exactly
    EXPECT_NEAR(rho_.get_residual(0, data_), 1.0, 1.0e-12);
}

TEST_F(DFPTRhoSerialTest, MixDrhoSecondStepCombinesCorrectly)
{
    psi::Psi<std::complex<double>> psi;
    std::vector<std::vector<std::vector<std::complex<double>>>> dpsi;
    FillRandomStates(psi, dpsi);
    ModuleBase::matrix wg(1, nbands_);
    wg(0, 0) = 1.0;
    wg(0, 1) = 0.0;

    rho_.compute_drho(psi, wg, 0, data_);
    rho_.mix_drho(0, data_);
    const std::vector<std::complex<double>> in1 = data_.get_drho_g(0, 0);

    // new response from fresh dpsi
    dfpt_rho_test_rng::g_seed = 777u;
    std::vector<std::vector<std::vector<std::complex<double>>>> dpsi2;
    FillRandomStates(psi, dpsi2);
    rho_.compute_drho(psi, wg, 0, data_);
    const std::vector<std::complex<double>> out2 = data_.get_drho_g(0, 0);
    rho_.mix_drho(0, data_);
    const std::vector<std::complex<double>> mixed2 = data_.get_drho_g(0, 0);

    double dn2 = 0.0;
    double o2 = 0.0;
    for (int ig = 0; ig < pw_rho_.npw; ++ig)
    {
        const std::complex<double> ref = in1[ig] + 0.4 * (out2[ig] - in1[ig]);
        EXPECT_NEAR(mixed2[ig].real(), ref.real(), 1.0e-12);
        EXPECT_NEAR(mixed2[ig].imag(), ref.imag(), 1.0e-12);
        dn2 += std::norm(out2[ig] - in1[ig]);
        o2 += std::norm(out2[ig]);
    }
    EXPECT_NEAR(rho_.get_residual(0, data_), std::sqrt(dn2 / o2), 1.0e-12);
}

TEST_F(DFPTRhoSerialTest, MixDrhoKerkerFirstStepIsPreconditionedScaledOutput)
{
    // first step from the zero input: mixed[ig] = beta * f[ig] * out[ig]
    // with the Kerker screen f[ig] = |G+q|^2 / (|G+q|^2 + a^2) built from
    // gcar + q_frac * G (1/lat0^2 units), the v_hartree_q convention
    ModuleDFPT::DFPT_Rho rho_k;
    double w2_min = 0.0;
    for (int ig = 0; ig < pw_rho_.npw; ++ig)
    {
        const ModuleBase::Vector3<double> w = pw_rho_.gcar[ig] + q_cart_;
        const double w2 = w * w;
        if (w2 > 1.0e-12 && (w2_min == 0.0 || w2 < w2_min))
        {
            w2_min = w2;
        }
    }
    ASSERT_GT(w2_min, 0.0);
    const double a2 = 4.0 * w2_min;
    rho_k.init({1, pw_rho_.nrxx, &pw_rho_, &pw_wfc_, G_, std::string("rho_kerker"), 0.7, a2});

    std::vector<std::complex<double>> out(pw_rho_.npw);
    int n_small = 0;
    int n_large = 0;
    for (int ig = 0; ig < pw_rho_.npw; ++ig)
    {
        out[ig] = 0.01 * std::complex<double>(std::cos(0.7 * ig), std::sin(0.5 * ig));
        const ModuleBase::Vector3<double> w = pw_rho_.gcar[ig] + q_cart_;
        const double w2 = w * w;
        if (w2 < a2)
        {
            ++n_small;
        }
        if (w2 > 100.0 * a2)
        {
            ++n_large;
        }
    }
    // sanity: the screen actually varies across the basis
    ASSERT_GT(n_small, 0);
    ASSERT_GT(n_large, 0);

    data_.set_drho_g(0, 0, out);
    rho_k.mix_drho(0, data_);
    const std::vector<std::complex<double>> mixed = data_.get_drho_g(0, 0);
    ASSERT_EQ(mixed.size(), out.size());
    for (int ig = 0; ig < pw_rho_.npw; ++ig)
    {
        const ModuleBase::Vector3<double> w = pw_rho_.gcar[ig] + q_cart_;
        const double w2 = w * w;
        const double f = (w2 < 1.0e-12) ? 0.0 : w2 / (w2 + a2);
        const std::complex<double> ref = 0.7 * f * out[ig];
        EXPECT_NEAR(mixed[ig].real(), ref.real(), 1.0e-12);
        EXPECT_NEAR(mixed[ig].imag(), ref.imag(), 1.0e-12);
    }
    EXPECT_NEAR(rho_k.get_residual(0, data_), 1.0, 1.0e-12);
}

TEST_F(DFPTRhoSerialTest, MixDrhoKerkerStabilizesStiffModelProblem)
{
    // model SCF problem: out(g) = D(g) in(g) + s(g) with the measured
    // diamond-smoke Coulomb-stiffness eigenvalue D = lambda ~ -2.2 on the
    // smallest |G+q| shell and D = 0.3 elsewhere; the fixed point is a
    // fixed target pattern t(g), s = (1 - D) t. Plain mixing must satisfy
    // beta < 2 / (1 + |lambda|) = 0.625, so beta = 0.7 diverges
    // (amplification |1 - beta (1 - D)| = 1.24), while the Kerker screen
    // damps the stiff shell amplification below 1 and converges.
    const int npw = pw_rho_.npw;

    double w2_min = 0.0;
    for (int ig = 0; ig < npw; ++ig)
    {
        const ModuleBase::Vector3<double> w = pw_rho_.gcar[ig] + q_cart_;
        const double w2 = w * w;
        if (w2 > 1.0e-12 && (w2_min == 0.0 || w2 < w2_min))
        {
            w2_min = w2;
        }
    }
    ASSERT_GT(w2_min, 0.0);

    std::vector<double> stiff(npw, 0.3);
    for (int ig = 0; ig < npw; ++ig)
    {
        const ModuleBase::Vector3<double> w = pw_rho_.gcar[ig] + q_cart_;
        const double w2 = w * w;
        if (w2 > 1.0e-12 && w2 < 1.5 * w2_min)
        {
            stiff[ig] = -2.2;
        }
    }

    std::vector<std::complex<double>> target(npw);
    for (int ig = 0; ig < npw; ++ig)
    {
        target[ig] = 0.01 * std::complex<double>(std::cos(0.3 * ig), std::sin(0.9 * ig));
    }

    auto model_out = [&](const std::vector<std::complex<double>>& in) {
        std::vector<std::complex<double>> o(npw);
        for (int ig = 0; ig < npw; ++ig)
        {
            o[ig] = stiff[ig] * in[ig] + (1.0 - stiff[ig]) * target[ig];
        }
        return o;
    };

    // plain beta = 0.7 on the stiff model diverges
    ModuleDFPT::DFPT_Rho rho_p;
    rho_p.init({1, pw_rho_.nrxx, &pw_rho_, &pw_wfc_, G_, std::string("rho_plain"), 0.7, 0.0});
    data_.set_drho_g(0, 0, std::vector<std::complex<double>>(npw, std::complex<double>(0.0, 0.0)));
    for (int it = 0; it < 40; ++it)
    {
        data_.set_drho_g(0, 0, model_out(data_.get_drho_g(0, 0)));
        rho_p.mix_drho(0, data_);
    }
    const double residual_plain = rho_p.get_residual(0, data_);
    EXPECT_GT(residual_plain, 1.0);

    // kerker beta = 0.7 with a^2 = 9 w2_min (f ~ 0.1 on the stiff shell)
    // converges to the target
    ModuleDFPT::DFPT_Rho rho_k;
    rho_k.init({1, pw_rho_.nrxx, &pw_rho_, &pw_wfc_, G_, std::string("rho_kerker"), 0.7, 9.0 * w2_min});
    data_.set_drho_g(0, 0, std::vector<std::complex<double>>(npw, std::complex<double>(0.0, 0.0)));
    for (int it = 0; it < 300; ++it)
    {
        data_.set_drho_g(0, 0, model_out(data_.get_drho_g(0, 0)));
        rho_k.mix_drho(0, data_);
    }
    const double residual_kerker = rho_k.get_residual(0, data_);
    EXPECT_LT(residual_kerker, 1.0e-8);
    const std::vector<std::complex<double>> final_in = data_.get_drho_g(0, 0);
    for (int ig = 0; ig < npw; ++ig)
    {
        EXPECT_NEAR(final_in[ig].real(), target[ig].real(), 1.0e-10);
        EXPECT_NEAR(final_in[ig].imag(), target[ig].imag(), 1.0e-10);
    }
}

TEST_F(DFPTRhoSerialTest, BroydenConvergesDiagonalSystemFast)
{
    // synthetic fixed point: out = D in + (1 - D) t with a diagonal D
    // (0.5 and 0.3 on two G shells, 0 elsewhere) and a fixed target t,
    // so the response map has a source and the relative residual
    // ||out - in|| / ||out|| genuinely decays. Plain beta = 0.4
    // contracts every populated shell by |1 - beta (1 - D)| = 0.8 (d = 0.5,
    // the binding shell; 0.72 for d = 0.3), i.e. roughly 100 plain
    // iterations to 1e-10; Modified Broyden on the few distinct
    // eigenvalues should converge in a handful of iterations.
    SetupBases(k_d_, q_d_, 2);
    ModuleDFPT::DFPT_Rho::Config cfg;
    cfg.nspin = 1;
    cfg.nrxx = pw_rho_.nrxx;
    cfg.pw_rho = &pw_rho_;
    cfg.pw_wfc = &pw_wfc_;
    cfg.recip_matrix = G_;
    cfg.mix_type = "rho_broyden";
    cfg.mix_beta = 0.4;
    cfg.kerker_a2 = 1.0;
    cfg.broyden_ndim = 8;
    ModuleDFPT::DFPT_Rho rho;
    rho.init(cfg);

    const int npw = pw_rho_.npw;
    std::vector<std::complex<double>> target(npw, std::complex<double>(0.0, 0.0));
    target[0] = std::complex<double>(0.9, -0.2);
    target[1] = std::complex<double>(0.5, 0.2);
    std::vector<double> d(npw, 0.0);
    d[0] = 0.5;
    d[1] = 0.3;
    auto model_out = [&](const std::vector<std::complex<double>>& in) {
        std::vector<std::complex<double>> o(npw);
        for (int ig = 0; ig < npw; ++ig)
        {
            o[ig] = d[ig] * in[ig] + (1.0 - d[ig]) * target[ig];
        }
        return o;
    };

    // run mixing iterations from the zero input; each step the model
    // rebuilds the output from the mixed input
    int n_iter = 0;
    std::vector<std::complex<double>> rho_out = model_out(std::vector<std::complex<double>>(npw, std::complex<double>(0.0, 0.0)));
    for (int iter = 0; iter < 40; ++iter)
    {
        data_.set_drho_g(0, 0, rho_out);
        rho.mix_drho(0, data_);
        ++n_iter;
        const double res = rho.get_residual(0, data_);
        if (res < 1.0e-10)
        {
            break;
        }
        rho_out = model_out(data_.get_drho_g(0, 0));
    }
    // Broyden on a 3-eigenvalue diagonal system: converges dramatically
    // faster than plain mixing (which needs roughly 100 iterations at
    // 0.8/step on the binding shell)
    EXPECT_LT(n_iter, 12);
    const double final_res = rho.get_residual(0, data_);
    EXPECT_LT(final_res, 1.0e-8);
    // the converged input is the fixed point of the model map
    const std::vector<std::complex<double>> final_in = data_.get_drho_g(0, 0);
    for (int ig = 0; ig < npw; ++ig)
    {
        EXPECT_NEAR(final_in[ig].real(), target[ig].real(), 1.0e-8);
        EXPECT_NEAR(final_in[ig].imag(), target[ig].imag(), 1.0e-8);
    }
}

TEST_F(DFPTRhoSerialTest, BroydenResetMixingAllowsFreshQAndSameQCycle)
{
    // regression (PR review): entering a fresh q after q = 0 must reset
    // the shared Broyden mixer's address binding even though the new q
    // has no history entry yet -- otherwise the first push of the newly
    // allocated Mixing_Data exits through the one-object-one-data guard
    // (broyden_mixing.cpp) and kills the process. The same-q cycle must
    // also restart from a zero input density.
    SetupBases(k_d_, q_d_, 2);
    ModuleDFPT::DFPT_Rho::Config cfg;
    cfg.nspin = 1;
    cfg.nrxx = pw_rho_.nrxx;
    cfg.pw_rho = &pw_rho_;
    cfg.pw_wfc = &pw_wfc_;
    cfg.recip_matrix = G_;
    cfg.mix_type = "rho_broyden";
    cfg.mix_beta = 0.4;
    cfg.kerker_a2 = 1.0;
    cfg.broyden_ndim = 8;
    ModuleDFPT::DFPT_Rho rho;
    rho.init(cfg);

    const int npw = pw_rho_.npw;
    std::vector<std::complex<double>> out(npw);
    for (int ig = 0; ig < npw; ++ig)
    {
        out[ig] = 0.01 * std::complex<double>(std::cos(0.3 * ig), std::sin(0.9 * ig));
    }

    // two pushes at q = 0: the mixer binds the q = 0 Mixing_Data
    data_.set_drho_g(0, 0, out);
    rho.mix_drho(0, data_);
    data_.set_drho_g(0, 0, out);
    rho.mix_drho(0, data_);

    // entering q = 1 for the first time (no history entry exists): with
    // the buggy guard the mixer stayed bound to q = 0 and this died
    rho.reset_mixing(1);
    data_.set_drho_g(1, 0, out);
    rho.mix_drho(1, data_);
    // fresh q starts from the zero input: first step is beta * out
    const std::vector<std::complex<double>> mixed_q1 = data_.get_drho_g(1, 0);
    for (int ig = 0; ig < npw; ++ig)
    {
        EXPECT_NEAR(mixed_q1[ig].real(), 0.4 * out[ig].real(), 1.0e-12);
        EXPECT_NEAR(mixed_q1[ig].imag(), 0.4 * out[ig].imag(), 1.0e-12);
    }

    // the same-q cycle also restarts from a zero input density
    rho.reset_mixing(0);
    data_.set_drho_g(0, 0, out);
    rho.mix_drho(0, data_);
    const std::vector<std::complex<double>> mixed_q0 = data_.get_drho_g(0, 0);
    for (int ig = 0; ig < npw; ++ig)
    {
        EXPECT_NEAR(mixed_q0[ig].real(), 0.4 * out[ig].real(), 1.0e-12);
        EXPECT_NEAR(mixed_q0[ig].imag(), 0.4 * out[ig].imag(), 1.0e-12);
    }
}
