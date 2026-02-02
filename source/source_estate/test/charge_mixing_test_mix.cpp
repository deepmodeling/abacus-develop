// Mix rho tests for charge_mixing
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#define private public
#include "../module_charge/charge_mixing.h"
#include "source_base/module_mixing/broyden_mixing.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_io/module_parameter/parameter.h"

class ChargeMixingMixTest : public ::testing::Test
{
  public:
    UnitCell ucell;
    ChargeMixingMixTest()
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

TEST_F(ChargeMixingMixTest, MixRhoTest)
{
     PARAM.sys.double_grid = false;
    charge.set_rhopw(&pw_basis);
    const int nspin = PARAM.input.nspin = 1;
    PARAM.sys.domag_z = false;
    XC_Functional::func_type = 3;
    XC_Functional::ked_flag = true;
    PARAM.input.mixing_beta = 0.7;
    PARAM.input.mixing_ndim = 1;
    PARAM.input.mixing_gg0 = 0.0;
    PARAM.input.mixing_tau = true;
    PARAM.input.mixing_mode = "plain";
    const int nrxx = pw_basis.nrxx;
    const int npw = pw_basis.npw;
    charge._space_rho = new double[nspin * nrxx];
    charge._space_rho_save = new double[nspin * nrxx];
    charge._space_rhog = new std::complex<double>[nspin * npw];
    charge._space_rhog_save = new std::complex<double>[nspin * npw];
    charge._space_kin_r = new double[nspin * nrxx];
    charge._space_kin_r_save = new double[nspin * nrxx];
    charge.rho = new double*[nspin];
    charge.rhog = new std::complex<double>*[nspin];
    charge.rho_save = new double*[nspin];
    charge.rhog_save = new std::complex<double>*[nspin];
    charge.kin_r = new double*[nspin];
    charge.kin_r_save = new double*[nspin];
    for (int is = 0; is < nspin; is++)
    {
        charge.rho[is] = charge._space_rho + is * nrxx;
        charge.rhog[is] = charge._space_rhog + is * npw;
        charge.rho_save[is] = charge._space_rho_save + is * nrxx;
        charge.rhog_save[is] = charge._space_rhog_save + is * npw;
        charge.kin_r[is] = charge._space_kin_r + is * nrxx;
        charge.kin_r_save[is] = charge._space_kin_r_save + is * nrxx;
    }
    std::vector<double> real_ref(nspin * nrxx);
    std::vector<double> real_save_ref(nspin * nrxx);
    std::vector<std::complex<double>> recip_ref(nspin * npw);
    std::vector<std::complex<double>> recip_save_ref(nspin * npw);
    for(int i = 0 ; i < nspin * npw; ++i)
    {
       recip_ref[i] = std::complex<double>(double(i), 1.0);
       recip_save_ref[i] = std::complex<double>(double(i), 0.0);
    }
    for(int i = 0 ; i < nspin ; ++i)
    {
        pw_basis.recip2real(recip_ref.data() + i * npw, real_ref.data() + i * nrxx);
        pw_basis.recip2real(recip_save_ref.data() + i * npw, real_save_ref.data() + i * nrxx);
    }
    //--------------------------------MAIN BODY--------------------------------
    // RECIPROCAL
    Charge_Mixing CMtest_recip;
    CMtest_recip.set_rhopw(&pw_basis, &pw_basis);
    PARAM.input.scf_thr_type= 1;
    CMtest_recip.set_mixing(PARAM.input.mixing_mode,
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
    CMtest_recip.init_mixing();
    for(int i = 0 ; i < nspin * npw; ++i)
    {
        charge._space_rhog[i] = recip_ref[i];
        charge._space_rhog_save[i] = recip_save_ref[i];
    }
    for(int i = 0 ; i < nspin * nrxx; ++i)
    {
        charge._space_rho[i] = real_ref[i];
        charge._space_rho_save[i] = real_save_ref[i];
    }
    CMtest_recip.mix_rho(&charge);
    for(int is = 0 ; is < nspin; ++is)
    {
        for(int ir = 0 ; ir < nrxx ; ++ir)
        {
            EXPECT_NEAR(charge.rho_save[is][ir], real_ref[is*nrxx + ir], 1e-8);
        }
        for(int ig = 0; ig < npw ; ++ig)
        {
            EXPECT_NEAR(charge.rhog[is][ig].real(), recip_save_ref[is*npw + ig].real(), 1e-8);
            EXPECT_NEAR(charge.rhog[is][ig].imag(), recip_save_ref[is*npw + ig].imag() + 0.7, 1e-8);
        }
    }

    // REAL
    Charge_Mixing CMtest_real;
    PARAM.input.scf_thr_type= 2;
    CMtest_real.set_rhopw(&pw_basis, &pw_basis);
    CMtest_real.set_mixing(PARAM.input.mixing_mode,
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
    CMtest_real.init_mixing();
    for(int i = 0 ; i < nspin * nrxx; ++i)
    {
        charge._space_rho[i] = real_ref[i];
        charge._space_rho_save[i] = real_save_ref[i];
    }
    CMtest_recip.mix_rho(&charge);
    for(int is = 0 ; is < nspin; ++is)
    {
        for(int ir = 0 ; ir < nrxx ; ++ir)
        {
            EXPECT_NEAR(charge.rho_save[is][ir], real_ref[is*nrxx + ir], 1e-8);
            EXPECT_NEAR(charge.rho[is][ir], 0.3*real_save_ref[is*nrxx+ir] + 0.7*real_ref[is*nrxx+ir], 1e-8);
        }
    }

    //-------------------------------------------------------------------------
    delete[] charge._space_rho;
    delete[] charge._space_rho_save;
    delete[] charge._space_rhog;
    delete[] charge._space_rhog_save;
    delete[] charge._space_kin_r;
    delete[] charge._space_kin_r_save;
    delete[] charge.rho;
    delete[] charge.rhog;
    delete[] charge.rho_save;
    delete[] charge.rhog_save;
    delete[] charge.kin_r;
    delete[] charge.kin_r_save;
}

TEST_F(ChargeMixingMixTest, MixDoubleGridRhoTest)
{
     PARAM.sys.double_grid = true;
    charge.set_rhopw(&pw_dbasis);
    const int nspin = PARAM.input.nspin = 1;
    PARAM.sys.domag_z = false;
    XC_Functional::func_type = 3;
    XC_Functional::ked_flag = true;
    PARAM.input.mixing_beta = 0.7;
    PARAM.input.mixing_ndim = 1;
    PARAM.input.mixing_gg0 = 0.0;
    PARAM.input.mixing_tau = true;
    PARAM.input.mixing_mode = "plain";
    const int nrxx = pw_dbasis.nrxx;
    const int npw = pw_dbasis.npw;
    charge._space_rho = new double[nspin * nrxx];
    charge._space_rho_save = new double[nspin * nrxx];
    charge._space_rhog = new std::complex<double>[nspin * npw];
    charge._space_rhog_save = new std::complex<double>[nspin * npw];
    charge._space_kin_r = new double[nspin * nrxx];
    charge._space_kin_r_save = new double[nspin * nrxx];
    charge.rho = new double*[nspin];
    charge.rhog = new std::complex<double>*[nspin];
    charge.rho_save = new double*[nspin];
    charge.rhog_save = new std::complex<double>*[nspin];
    charge.kin_r = new double*[nspin];
    charge.kin_r_save = new double*[nspin];
    for (int is = 0; is < nspin; is++)
    {
        charge.rho[is] = charge._space_rho + is * nrxx;
        charge.rhog[is] = charge._space_rhog + is * npw;
        charge.rho_save[is] = charge._space_rho_save + is * nrxx;
        charge.rhog_save[is] = charge._space_rhog_save + is * npw;
        charge.kin_r[is] = charge._space_kin_r + is * nrxx;
        charge.kin_r_save[is] = charge._space_kin_r_save + is * nrxx;
    }
    std::vector<double> real_ref(nspin * nrxx);
    std::vector<double> real_save_ref(nspin * nrxx);
    std::vector<std::complex<double>> recip_ref(nspin * npw);
    std::vector<std::complex<double>> recip_save_ref(nspin * npw);
    for (int i = 0; i < nspin * npw; ++i)
    {
        recip_ref[i] = std::complex<double>(double(i), 1.0);
        recip_save_ref[i] = std::complex<double>(double(i), 0.0);
    }
    for (int i = 0; i < nspin; ++i)
    {
        pw_dbasis.recip2real(recip_ref.data() + i * npw, real_ref.data() + i * nrxx);
        pw_dbasis.recip2real(recip_save_ref.data() + i * npw, real_save_ref.data() + i * nrxx);
    }
    //--------------------------------MAIN BODY--------------------------------
    // RECIPROCAL
    Charge_Mixing CMtest_recip;
    CMtest_recip.set_rhopw(&pw_basis, &pw_dbasis);

    PARAM.input.scf_thr_type= 1;
    CMtest_recip.set_mixing(PARAM.input.mixing_mode,
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

    CMtest_recip.init_mixing();
    for (int i = 0; i < nspin * npw; ++i)
    {
        charge._space_rhog[i] = recip_ref[i];
        charge._space_rhog_save[i] = recip_save_ref[i];
    }
    for (int i = 0; i < nspin * nrxx; ++i)
    {
        charge._space_rho[i] = real_ref[i];
        charge._space_rho_save[i] = real_save_ref[i];
    }
    CMtest_recip.mix_rho(&charge);
    for (int is = 0; is < nspin; ++is)
    {
        for (int ir = 0; ir < nrxx; ++ir)
        {
            EXPECT_NEAR(charge.rho_save[is][ir], real_ref[is * nrxx + ir], 1e-8);
        }
        for (int ig = 0; ig < npw; ++ig)
        {
            EXPECT_NEAR(charge.rhog[is][ig].real(), recip_save_ref[is * npw + ig].real(), 1e-8);
            EXPECT_NEAR(charge.rhog[is][ig].imag(), recip_save_ref[is * npw + ig].imag() + 0.7, 1e-8);
        }
    }

    //-------------------------------------------------------------------------
    delete[] charge._space_rho;
    delete[] charge._space_rho_save;
    delete[] charge._space_rhog;
    delete[] charge._space_rhog_save;
    delete[] charge._space_kin_r;
    delete[] charge._space_kin_r_save;
    delete[] charge.rho;
    delete[] charge.rhog;
    delete[] charge.rho_save;
    delete[] charge.rhog_save;
    delete[] charge.kin_r;
    delete[] charge.kin_r_save;
}

TEST_F(ChargeMixingMixTest, MixDivCombTest)
{
    // NSPIN = 1
    PARAM.input.nspin = 1;
    Charge_Mixing CMtest;
    CMtest.set_rhopw(&pw_basis, &pw_dbasis);
    std::vector<std::complex<double>> data(pw_dbasis.npw, 1.0);
    std::complex<double>*datas, *datahf;
    std::complex<double>*datas2, *datahf2;
    CMtest.divide_data(data.data(), datas, datahf);
    EXPECT_EQ(datas, data.data());
    EXPECT_EQ(datahf, data.data() + pw_basis.npw);
    CMtest.combine_data(data.data(), datas, datahf);
    EXPECT_EQ(datas, nullptr);
    EXPECT_EQ(datahf, nullptr);

    CMtest.divide_data(data.data(), datas2, datahf2);
    CMtest.clean_data(datas2, datahf2);
    EXPECT_EQ(datas2, nullptr);
    EXPECT_EQ(datahf2, nullptr);

    // NSPIN = 2
    PARAM.input.nspin = 2;
    data.resize(pw_dbasis.npw * 2, 1.0);
    std::vector<std::complex<double>> dataout(pw_dbasis.npw * 2, 1.0);
    CMtest.divide_data(data.data(), datas, datahf);
    CMtest.combine_data(dataout.data(), datas, datahf);
    EXPECT_EQ(datas, nullptr);
    EXPECT_EQ(datahf, nullptr);
    for (int i = 0; i < pw_dbasis.npw * 2; ++i)
    {
        EXPECT_EQ(dataout[i], data[i]);
    }

    CMtest.divide_data(data.data(), datas2, datahf2);
    CMtest.clean_data(datas2, datahf2);
    EXPECT_EQ(datas2, nullptr);
    EXPECT_EQ(datahf2, nullptr);
}

TEST_F(ChargeMixingMixTest, SCFOscillationTest)
{
    Charge_Mixing CMtest;
    int scf_nmax = 20;
    int scf_os_ndim = 3;
    double scf_os_thr = -0.05;
    bool scf_oscillate = false;
    std::vector<double> drho(scf_nmax, 0.0);
    std::vector<bool> scf_oscillate_ref(scf_nmax, false);
    drho = {6.83639633652e-05,
            4.93523029235e-05,
            3.59230097735e-05,
            2.68356403913e-05,
            2.17490806464e-05,
            2.14231642508e-05,
            1.67507494811e-05,
            1.53575889539e-05,
            1.26504511554e-05,
            1.04762016224e-05,
            8.10000162918e-06,
            7.66427917682e-06,
            6.70112820094e-06,
            5.68594436664e-06,
            4.80120233733e-06,
            4.86519757184e-06,
            4.37855804356e-06,
            4.29922703412e-06,
            4.36398486331e-06,
            4.94224615955e-06};
    scf_oscillate_ref = {false,false,false,false,false,true,false,false,false,false,
                        false,false,true,false,false,true,true,true,true,true};
    for (int i = 1; i <= scf_nmax; ++i)
    {
        scf_oscillate = CMtest.if_scf_oscillate(i,drho[i-1],scf_os_ndim,scf_os_thr);
        EXPECT_EQ(scf_oscillate, scf_oscillate_ref[i-1]);
    }
}
