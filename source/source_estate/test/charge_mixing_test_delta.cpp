// Delta rho mixing tests for charge_mixing (DFPT)
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#define private public
#include "../module_charge/charge_mixing.h"
#include "source_base/module_mixing/broyden_mixing.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_io/module_parameter/parameter.h"

class ChargeMixingDeltaTest : public ::testing::Test
{
  public:
    UnitCell ucell;
    ChargeMixingDeltaTest()
    {
        // Init pw_basis
        pw_basis.initgrids(4, ModuleBase::Matrix3(1, 0, 0, 0, 1, 0, 0, 0, 1), 20);
        pw_basis.initparameters(false, 20);
        pw_basis.setuptransform();
        pw_basis.collect_local_pw();
        pw_dbasis.initgrids(4, ModuleBase::Matrix3(1, 0, 0, 0, 1, 0, 0, 0, 1), 40);
        pw_dbasis.initparameters(false, 40);
        pw_dbasis.setuptransform(&pw_basis);
        pw_dbasis.collect_local_pw();
        // default mixing parameters
        PARAM.input.mixing_mode = "broyden";
        PARAM.input.mixing_beta = 0.8;
        PARAM.input.mixing_ndim = 8;
        PARAM.input.mixing_gg0  = 1.0;
        PARAM.input.mixing_tau  = false;
        PARAM.input.mixing_beta_mag = 1.6;
        PARAM.input.mixing_gg0_mag = 0.0;
        PARAM.input.mixing_gg0_min = 0.1;
        PARAM.input.mixing_angle = -10.0;
        PARAM.input.mixing_dmr = false;
        ucell.omega = 1.0;
        ucell.tpiba = 1.0;
    }
    ModulePW::PW_Basis pw_basis;
    ModulePW::PW_Basis_Sup pw_dbasis;
    Charge charge;
};

TEST_F(ChargeMixingDeltaTest, MixDeltaRhoBroydenTest)
{
    // Test that delta_rho can be mixed correctly with Broyden mixing method
    // This is used for DFPT calculations where we mix complex charge density response

    const int nspin = 1;
    PARAM.input.nspin = nspin;
    const int nrxx = pw_basis.nrxx;
    const int npw = pw_basis.npw;

    // Setup mixing parameters for Broyden mixing
    PARAM.input.mixing_mode = "broyden";
    PARAM.input.mixing_beta = 0.7;
    PARAM.input.mixing_ndim = 8;
    PARAM.input.mixing_gg0 = 0.0;  // No Kerker screening for simplicity

    Charge_Mixing CMtest;
    CMtest.set_rhopw(&pw_basis, &pw_basis);
    CMtest.set_mixing(PARAM.input.mixing_mode,
                      PARAM.input.mixing_beta,
                      PARAM.input.mixing_ndim,
                      PARAM.input.mixing_gg0,
                      PARAM.input.mixing_tau,
                      PARAM.input.mixing_beta_mag,
                      PARAM.input.mixing_gg0_mag,
                      PARAM.input.mixing_gg0_min,
                      PARAM.input.mixing_angle,
                      PARAM.input.mixing_dmr,
                      ucell.omega,
                      ucell.tpiba);
    CMtest.init_mixing();

    // Initialize delta_rho mixing data
    CMtest.init_mixing_delta_rho(nspin * npw);

    // Allocate delta_rho arrays in Charge object
    charge._space_delta_rho = new std::complex<double>[nspin * nrxx];
    charge._space_delta_rho_save = new std::complex<double>[nspin * nrxx];
    charge._space_delta_rhog = new std::complex<double>[nspin * npw];
    charge._space_delta_rhog_save = new std::complex<double>[nspin * npw];
    charge.delta_rho = new std::complex<double>*[nspin];
    charge.delta_rho_save = new std::complex<double>*[nspin];
    charge.delta_rhog = new std::complex<double>*[nspin];
    charge.delta_rhog_save = new std::complex<double>*[nspin];
    for (int is = 0; is < nspin; is++)
    {
        charge.delta_rho[is] = charge._space_delta_rho + is * nrxx;
        charge.delta_rho_save[is] = charge._space_delta_rho_save + is * nrxx;
        charge.delta_rhog[is] = charge._space_delta_rhog + is * npw;
        charge.delta_rhog_save[is] = charge._space_delta_rhog_save + is * npw;
    }

    // Initialize test data: delta_rho (current) and delta_rho_save (previous)
    // Use simple patterns that can be verified after mixing
    std::vector<std::complex<double>> delta_rhog_in(nspin * npw);
    std::vector<std::complex<double>> delta_rhog_out(nspin * npw);
    for (int i = 0; i < nspin * npw; ++i)
    {
        // Previous iteration (input to mixing)
        delta_rhog_in[i] = std::complex<double>(1.0, 0.0);
        // Current iteration (output from SCF, to be mixed)
        delta_rhog_out[i] = std::complex<double>(2.0, 0.5);
    }

    // Convert to real space for the Charge object
    for (int is = 0; is < nspin; ++is)
    {
        pw_basis.recip2real(delta_rhog_in.data() + is * npw, charge.delta_rho_save[is]);
        pw_basis.recip2real(delta_rhog_out.data() + is * npw, charge.delta_rho[is]);
        for (int ig = 0; ig < npw; ++ig)
        {
            charge.delta_rhog_save[is][ig] = delta_rhog_in[is * npw + ig];
            charge.delta_rhog[is][ig] = delta_rhog_out[is * npw + ig];
        }
    }

    // Store original values for comparison
    std::vector<std::complex<double>> delta_rho_before(nspin * nrxx);
    for (int is = 0; is < nspin; ++is)
    {
        for (int ir = 0; ir < nrxx; ++ir)
        {
            delta_rho_before[is * nrxx + ir] = charge.delta_rho[is][ir];
        }
    }

    // Perform mixing - first iteration with plain mixing (Broyden needs history)
    CMtest.mix_delta_rho(&charge, nspin);

    // Verify that mixing occurred:
    // For the first iteration of Broyden, it should behave like plain mixing:
    // mixed = in + beta * (out - in) = in + beta * residual
    // In reciprocal space: mixed_g = delta_rhog_in + 0.7 * (delta_rhog_out - delta_rhog_in)
    //                              = 1.0 + 0.7 * (2.0 + 0.5i - 1.0)
    //                              = 1.0 + 0.7 * (1.0 + 0.5i)
    //                              = 1.7 + 0.35i
    for (int is = 0; is < nspin; ++is)
    {
        for (int ig = 0; ig < npw; ++ig)
        {
            std::complex<double> expected(1.7, 0.35);
            EXPECT_NEAR(charge.delta_rhog[is][ig].real(), expected.real(), 1e-8);
            EXPECT_NEAR(charge.delta_rhog[is][ig].imag(), expected.imag(), 1e-8);
        }
    }

    // Verify delta_rho_save was updated with pre-mixing values
    for (int is = 0; is < nspin; ++is)
    {
        for (int ir = 0; ir < nrxx; ++ir)
        {
            EXPECT_NEAR(charge.delta_rho_save[is][ir].real(),
                        delta_rho_before[is * nrxx + ir].real(), 1e-8);
            EXPECT_NEAR(charge.delta_rho_save[is][ir].imag(),
                        delta_rho_before[is * nrxx + ir].imag(), 1e-8);
        }
    }

    // Cleanup
    delete[] charge._space_delta_rho;
    delete[] charge._space_delta_rho_save;
    delete[] charge._space_delta_rhog;
    delete[] charge._space_delta_rhog_save;
    delete[] charge.delta_rho;
    delete[] charge.delta_rho_save;
    delete[] charge.delta_rhog;
    delete[] charge.delta_rhog_save;
    charge._space_delta_rho = nullptr;
    charge._space_delta_rho_save = nullptr;
    charge._space_delta_rhog = nullptr;
    charge._space_delta_rhog_save = nullptr;
    charge.delta_rho = nullptr;
    charge.delta_rho_save = nullptr;
    charge.delta_rhog = nullptr;
    charge.delta_rhog_save = nullptr;
}

TEST_F(ChargeMixingDeltaTest, MixDeltaRhoBroydenMultiIterTest)
{
    // Test Broyden mixing with multiple iterations to verify history-based mixing

    const int nspin = 1;
    PARAM.input.nspin = nspin;
    const int nrxx = pw_basis.nrxx;
    const int npw = pw_basis.npw;

    // Setup mixing parameters for Broyden mixing
    PARAM.input.mixing_mode = "broyden";
    PARAM.input.mixing_beta = 0.7;
    PARAM.input.mixing_ndim = 4;
    PARAM.input.mixing_gg0 = 0.0;

    Charge_Mixing CMtest;
    CMtest.set_rhopw(&pw_basis, &pw_basis);
    CMtest.set_mixing(PARAM.input.mixing_mode,
                      PARAM.input.mixing_beta,
                      PARAM.input.mixing_ndim,
                      PARAM.input.mixing_gg0,
                      PARAM.input.mixing_tau,
                      PARAM.input.mixing_beta_mag,
                      PARAM.input.mixing_gg0_mag,
                      PARAM.input.mixing_gg0_min,
                      PARAM.input.mixing_angle,
                      PARAM.input.mixing_dmr,
                      ucell.omega,
                      ucell.tpiba);
    CMtest.init_mixing();
    CMtest.init_mixing_delta_rho(nspin * npw);

    // Allocate delta_rho arrays
    charge._space_delta_rho = new std::complex<double>[nspin * nrxx];
    charge._space_delta_rho_save = new std::complex<double>[nspin * nrxx];
    charge._space_delta_rhog = new std::complex<double>[nspin * npw];
    charge._space_delta_rhog_save = new std::complex<double>[nspin * npw];
    charge.delta_rho = new std::complex<double>*[nspin];
    charge.delta_rho_save = new std::complex<double>*[nspin];
    charge.delta_rhog = new std::complex<double>*[nspin];
    charge.delta_rhog_save = new std::complex<double>*[nspin];
    for (int is = 0; is < nspin; is++)
    {
        charge.delta_rho[is] = charge._space_delta_rho + is * nrxx;
        charge.delta_rho_save[is] = charge._space_delta_rho_save + is * nrxx;
        charge.delta_rhog[is] = charge._space_delta_rhog + is * npw;
        charge.delta_rhog_save[is] = charge._space_delta_rhog_save + is * npw;
    }

    // Simulate multiple SCF iterations with decreasing residuals
    // This tests that Broyden mixing properly uses history
    std::vector<double> residuals;

    for (int iter = 0; iter < 3; ++iter)
    {
        // Set up input (previous) and output (current) for this iteration
        double scale_in = 1.0 + 0.5 * iter;
        double scale_out = 2.0 + 0.3 * iter;

        for (int is = 0; is < nspin; ++is)
        {
            for (int ig = 0; ig < npw; ++ig)
            {
                charge.delta_rhog_save[is][ig] = std::complex<double>(scale_in, 0.1 * iter);
                charge.delta_rhog[is][ig] = std::complex<double>(scale_out, 0.2 * iter);
            }
            pw_basis.recip2real(charge.delta_rhog_save[is], charge.delta_rho_save[is]);
            pw_basis.recip2real(charge.delta_rhog[is], charge.delta_rho[is]);
        }

        // Get residual before mixing
        double drho = CMtest.get_delta_drho(&charge, nspin);
        residuals.push_back(drho);

        // Perform mixing
        CMtest.mix_delta_rho(&charge, nspin);
    }

    // Verify that residuals were computed (non-zero for non-converged case)
    for (int i = 0; i < 3; ++i)
    {
        EXPECT_GT(residuals[i], 0.0);
    }

    // Verify that mixing data structure is properly maintained
    EXPECT_NE(CMtest.delta_rho_mdata.data, nullptr);
    EXPECT_EQ(CMtest.delta_rho_mdata.length, nspin * npw);

    // Cleanup
    delete[] charge._space_delta_rho;
    delete[] charge._space_delta_rho_save;
    delete[] charge._space_delta_rhog;
    delete[] charge._space_delta_rhog_save;
    delete[] charge.delta_rho;
    delete[] charge.delta_rho_save;
    delete[] charge.delta_rhog;
    delete[] charge.delta_rhog_save;
    charge._space_delta_rho = nullptr;
    charge._space_delta_rho_save = nullptr;
    charge._space_delta_rhog = nullptr;
    charge._space_delta_rhog_save = nullptr;
    charge.delta_rho = nullptr;
    charge.delta_rho_save = nullptr;
    charge.delta_rhog = nullptr;
    charge.delta_rhog_save = nullptr;
}

TEST_F(ChargeMixingDeltaTest, InnerProductRecipComplexTest)
{
    // Test the complex inner product function used for DFPT mixing

    Charge_Mixing CMtest;
    CMtest.set_rhopw(&pw_basis, &pw_basis);
    PARAM.input.nspin = 1;

    const int npw = pw_basis.npw;

    // Test 1: Inner product of vector with itself should give sum of |z|^2
    std::vector<std::complex<double>> rho1(npw);
    for (int i = 0; i < npw; ++i)
    {
        rho1[i] = std::complex<double>(1.0, 1.0);
    }
    double inner = CMtest.inner_product_recip_complex(rho1.data(), rho1.data());
    // |1+i|^2 = 2, so sum = 2 * npw
    EXPECT_NEAR(inner, 2.0 * npw, 1e-8);

    // Test 2: Inner product with orthogonal vectors
    std::vector<std::complex<double>> rho2(npw);
    for (int i = 0; i < npw; ++i)
    {
        rho1[i] = std::complex<double>(1.0, 0.0);
        rho2[i] = std::complex<double>(0.0, 1.0);
    }
    inner = CMtest.inner_product_recip_complex(rho1.data(), rho2.data());
    // conj(1) * i = 1 * i = i, real part = 0
    EXPECT_NEAR(inner, 0.0, 1e-8);

    // Test 3: General case
    for (int i = 0; i < npw; ++i)
    {
        rho1[i] = std::complex<double>(2.0, 1.0);
        rho2[i] = std::complex<double>(1.0, 2.0);
    }
    inner = CMtest.inner_product_recip_complex(rho1.data(), rho2.data());
    // conj(2+i) * (1+2i) = (2-i) * (1+2i) = 2 + 4i - i - 2i^2 = 2 + 3i + 2 = 4 + 3i
    // real part = 4, sum = 4 * npw
    EXPECT_NEAR(inner, 4.0 * npw, 1e-8);
}
