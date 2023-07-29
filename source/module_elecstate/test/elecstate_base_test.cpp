#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#define private public
#define protected public
#include "module_elecstate/elecstate.h"
#include "module_elecstate/occupy.h"

// Mock functions for testing elecstate.cpp
namespace elecstate
{
void Potential::init_pot(int, Charge const*)
{
}
void Potential::cal_v_eff(const Charge* chg, const UnitCell* ucell, ModuleBase::matrix& v_eff)
{
}
void Potential::cal_fixed_v(double* vl_pseudo)
{
}
Potential::~Potential()
{
}
} // namespace elecstate
Charge::Charge()
{
}
Charge::~Charge()
{
}
K_Vectors::K_Vectors()
{
}
K_Vectors::~K_Vectors()
{
}
ModulePW::PW_Basis::PW_Basis()
{
}
ModulePW::PW_Basis::~PW_Basis()
{
}
ModulePW::FFT::FFT()
{
}
ModulePW::FFT::~FFT()
{
}
void ModulePW::PW_Basis::initgrids(double, ModuleBase::Matrix3, double)
{
}
void ModulePW::PW_Basis::initgrids(double, ModuleBase::Matrix3, int, int, int)
{
}
void ModulePW::PW_Basis::distribute_r()
{
}
void Charge::set_rho_core(ModuleBase::ComplexMatrix const&)
{
}
void Charge::init_rho(elecstate::efermi&, ModuleBase::ComplexMatrix const&, const int&, const int&)
{
}
void Charge::set_rhopw(ModulePW::PW_Basis*)
{
}
void Charge::renormalize_rho()
{
}

/************************************************
 *  unit test of elecstate.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - InitNelecSpin: elecstate::ElecState::init_nelec_spin()
 *      - determine the number of electrons for spin up and down
 *   - Constructor: elecstate::ElecState(charge, rhopw, bigpw)
 *      - constructor ElecState using existing charge, rhopw, bigpw
 *   - CalNbands: elecstate::ElecState::cal_nbands()
 *      - calculate the number of bands
 *   - InitKS: elecstate::ElecState::init_ks()
 *      - initialize the elecstate for KS-DFT
 *   - GetRho: elecstate::ElecState::getRho()
 *      - get the pointer to this->charge->rho
 *   - VirtualBaseFuncs:
 *      - trivial calling to virtual functions including
 *         - elecstate::ElecState::psiToRho()
 *         - elecstate::ElecState::print_psi()
 *         - elecstate::ElecState::getNewRho()
 *    - InitSCF: elecstate::ElecState::init_scf()
 *      - trivial calling to elecstate::ElecState::init_scf()
 *      - the production function is init charge and pot for scf calculation
 *    - FixedWeights: elecstate::ElecState::fixed_weights()
 *      - fix wg using external weights: ocp_kb
 *    - CalEBand: elecstate::ElecState::cal_eband()
 *      - calculate the electronic states energy contribution to total energy
 *    - CalculateWeights: elecstate::ElecState::calculate_weights()
 *     - calculate the weights for each electronic state
 */

namespace elecstate
{
class MockElecState : public ElecState
{
  public:
    void Set_GlobalV_Default()
    {
        GlobalV::NSPIN = 1;
        GlobalV::nelec = 10.0;
        GlobalV::nupdown = 0.0;
        GlobalV::TWO_EFERMI = false;
        GlobalV::NBANDS = 6;
        GlobalV::NLOCAL = 6;
        GlobalV::ESOLVER_TYPE = "ksdft";
        GlobalV::LSPINORB = false;
        GlobalV::BASIS_TYPE = "pw";
        GlobalV::KPAR = 1;
        GlobalV::NPROC_IN_POOL = 1;
    }
};
} // namespace elecstate

class ElecStateTest : public ::testing::Test
{
  protected:
    elecstate::MockElecState* elecstate;
    std::string output;
    void SetUp()
    {
        elecstate = new elecstate::MockElecState;
        elecstate->Set_GlobalV_Default();
    }
    void TearDown()
    {
        delete elecstate;
    }
};

using ElecStateDeathTest = ElecStateTest;

TEST_F(ElecStateTest, InitNelecSpin)
{
    GlobalV::NSPIN = 2;
    elecstate->init_nelec_spin();
    EXPECT_EQ(elecstate->nelec_spin[0], 5.0);
    EXPECT_EQ(elecstate->nelec_spin[1], 5.0);
}

TEST_F(ElecStateTest, Constructor)
{
    Charge* charge = new Charge;
    ModulePW::PW_Basis* rhopw = new ModulePW::PW_Basis;
    ModulePW::PW_Basis_Big* bigpw = new ModulePW::PW_Basis_Big;
    elecstate::ElecState* elecstate_new = new elecstate::ElecState(charge, rhopw, bigpw);
    EXPECT_EQ(elecstate_new->charge, charge);
    EXPECT_EQ(elecstate_new->bigpw, bigpw);
    EXPECT_EQ(elecstate_new->eferm.two_efermi, GlobalV::TWO_EFERMI);
    delete elecstate_new;
    delete bigpw;
    delete rhopw;
    delete charge;
}

TEST_F(ElecStateTest, CalNbands)
{
    elecstate->cal_nbands();
    EXPECT_EQ(GlobalV::NBANDS, 6);
}

TEST_F(ElecStateTest, CalNbandsFractionElec)
{
    GlobalV::nelec = 9.5;
    elecstate->cal_nbands();
    EXPECT_EQ(GlobalV::NBANDS, 6);
}

TEST_F(ElecStateTest, CalNbandsSOC)
{
    GlobalV::LSPINORB = true;
    GlobalV::NBANDS = 0;
    elecstate->cal_nbands();
    EXPECT_EQ(GlobalV::NBANDS, 20);
}

TEST_F(ElecStateTest, CalNbandsSDFT)
{
    GlobalV::ESOLVER_TYPE = "sdft";
    EXPECT_NO_THROW(elecstate->cal_nbands());
}

TEST_F(ElecStateTest, CalNbandsLCAO)
{
    GlobalV::BASIS_TYPE = "lcao";
    EXPECT_NO_THROW(elecstate->cal_nbands());
}

TEST_F(ElecStateDeathTest, CalNbandsLCAOINPW)
{
    GlobalV::BASIS_TYPE = "lcao_in_pw";
    GlobalV::NLOCAL = GlobalV::NBANDS - 1;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(elecstate->cal_nbands(), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("NLOCAL < NBANDS"));
}

TEST_F(ElecStateDeathTest, CalNbandsWarning1)
{
    GlobalV::NBANDS = GlobalV::nelec / 2 - 1;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(elecstate->cal_nbands(), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("Too few bands!"));
}

TEST_F(ElecStateDeathTest, CalNbandsWarning2)
{
    GlobalV::NSPIN = 2;
    GlobalV::nupdown = 4.0;
    elecstate->init_nelec_spin();
    testing::internal::CaptureStdout();
    EXPECT_EXIT(elecstate->cal_nbands(), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("Too few spin up bands!"));
}

TEST_F(ElecStateDeathTest, CalNbandsWarning3)
{
    GlobalV::NSPIN = 2;
    GlobalV::nupdown = -4.0;
    elecstate->init_nelec_spin();
    testing::internal::CaptureStdout();
    EXPECT_EXIT(elecstate->cal_nbands(), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("Too few spin down bands!"));
}

TEST_F(ElecStateTest, CalNbandsSpin1)
{
    GlobalV::NSPIN = 1;
    GlobalV::NBANDS = 0;
    elecstate->cal_nbands();
    EXPECT_EQ(GlobalV::NBANDS, 15);
}

TEST_F(ElecStateTest, CalNbandsSpin1LCAO)
{
    GlobalV::NSPIN = 1;
    GlobalV::NBANDS = 0;
    GlobalV::BASIS_TYPE = "lcao";
    elecstate->cal_nbands();
    EXPECT_EQ(GlobalV::NBANDS, 6);
}

TEST_F(ElecStateTest, CalNbandsSpin4)
{
    GlobalV::NSPIN = 4;
    GlobalV::NBANDS = 0;
    elecstate->cal_nbands();
    EXPECT_EQ(GlobalV::NBANDS, 30);
}

TEST_F(ElecStateTest, CalNbandsSpin4LCAO)
{
    GlobalV::NSPIN = 4;
    GlobalV::NBANDS = 0;
    GlobalV::BASIS_TYPE = "lcao";
    elecstate->cal_nbands();
    EXPECT_EQ(GlobalV::NBANDS, 6);
}

TEST_F(ElecStateTest, CalNbandsSpin2)
{
    GlobalV::NSPIN = 2;
    GlobalV::NBANDS = 0;
    elecstate->init_nelec_spin();
    elecstate->cal_nbands();
    EXPECT_EQ(GlobalV::NBANDS, 16);
}

TEST_F(ElecStateTest, CalNbandsSpin2LCAO)
{
    GlobalV::NSPIN = 2;
    GlobalV::NBANDS = 0;
    GlobalV::BASIS_TYPE = "lcao";
    elecstate->init_nelec_spin();
    elecstate->cal_nbands();
    EXPECT_EQ(GlobalV::NBANDS, 6);
}

TEST_F(ElecStateDeathTest, CalNbandsGaussWarning)
{
    Occupy::use_gaussian_broadening = true;
    EXPECT_TRUE(Occupy::gauss());
    GlobalV::NBANDS = 5;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(elecstate->cal_nbands(), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("for smearing, num. of bands > num. of occupied bands"));
    Occupy::use_gaussian_broadening = false;
}

TEST_F(ElecStateTest, InitKS)
{
    Charge* charge = new Charge;
    ModulePW::PW_Basis* rhopw = new ModulePW::PW_Basis;
    ModulePW::PW_Basis_Big* bigpw = new ModulePW::PW_Basis_Big;
    K_Vectors* klist = new K_Vectors;
    int nk = 1;
    EXPECT_NO_THROW(elecstate->init_ks(charge, klist, nk, rhopw, bigpw));
    EXPECT_EQ(elecstate->charge, charge);
    EXPECT_EQ(elecstate->bigpw, bigpw);
    EXPECT_EQ(elecstate->klist, klist);
    EXPECT_EQ(elecstate->ekb.nr, nk);
    EXPECT_EQ(elecstate->ekb.nc, GlobalV::NBANDS);
    EXPECT_EQ(elecstate->wg.nr, nk);
    EXPECT_EQ(elecstate->wg.nc, GlobalV::NBANDS);
    delete klist;
    delete bigpw;
    delete rhopw;
    delete charge;
}

TEST_F(ElecStateTest, GetRho)
{
    Charge* charge = new Charge;
    ModulePW::PW_Basis* rhopw = new ModulePW::PW_Basis;
    ModulePW::PW_Basis_Big* bigpw = new ModulePW::PW_Basis_Big;
    K_Vectors* klist = new K_Vectors;
    int nk = 1;
    int nrxx = 100;
    charge->rho = new double*[GlobalV::NSPIN];
    for (int i = 0; i < GlobalV::NSPIN; ++i)
    {
        charge->rho[i] = new double[nrxx];
        for (int j = 0; j < nrxx; ++j)
        {
            charge->rho[i][j] = 1.0;
        }
    }
    elecstate->init_ks(charge, klist, nk, rhopw, bigpw);
    EXPECT_EQ(elecstate->getRho(0), &(charge->rho[0][0]));
    EXPECT_EQ(elecstate->getRho(0)[nrxx - 1], 1.0);
    for (int i = 0; i < GlobalV::NSPIN; ++i)
    {
        delete[] charge->rho[i];
    }
    delete[] charge->rho;
    delete klist;
    delete bigpw;
    delete rhopw;
    delete charge;
}

TEST_F(ElecStateTest, VirtualBaseFuncs)
{
    psi::Psi<std::complex<double>> psi_complex;
    psi::Psi<double> psi_real;
    EXPECT_NO_THROW(elecstate->psiToRho(psi_complex));
    EXPECT_NO_THROW(elecstate->psiToRho(psi_real));
    EXPECT_NO_THROW(elecstate->print_psi(psi_complex));
    EXPECT_NO_THROW(elecstate->print_psi(psi_real));
    EXPECT_NO_THROW(elecstate->getNewRho());
}

TEST_F(ElecStateTest, InitSCF)
{
    Charge* charge = new Charge;
    elecstate->charge = charge;
    elecstate->pot = new elecstate::Potential;
    elecstate::efermi efermi;
    int istep = 0;
    ModuleBase::ComplexMatrix strucfac;
    elecstate->eferm = efermi;
    EXPECT_NO_THROW(elecstate->init_scf(istep, strucfac));
    // delete elecstate->pot is done in the destructor of elecstate
    delete charge;
}

TEST_F(ElecStateTest,FixedWeights)
{
    EXPECT_EQ(GlobalV::NBANDS, 6);
    GlobalV::nelec = 30;
    K_Vectors* klist = new K_Vectors;
    klist->nks = 5;
    elecstate->klist = klist;
    elecstate->wg.create(klist->nks, GlobalV::NBANDS);
    std::vector<double> ocp_kb;
    ocp_kb.resize(GlobalV::NBANDS*elecstate->klist->nks);
    for (int i = 0; i < ocp_kb.size(); ++i)
    {
        ocp_kb[i] = 1.0;
    }
    elecstate->fixed_weights(ocp_kb);
    EXPECT_EQ(elecstate->wg(0, 0), 1.0);
    EXPECT_EQ(elecstate->wg(klist->nks-1, GlobalV::NBANDS-1), 1.0);
    EXPECT_TRUE(elecstate->skip_weights);
}

TEST_F(ElecStateDeathTest,FixedWeightsWarning1)
{
    EXPECT_EQ(GlobalV::NBANDS, 6);
    GlobalV::nelec = 30;
    K_Vectors* klist = new K_Vectors;
    klist->nks = 5;
    elecstate->klist = klist;
    elecstate->wg.create(klist->nks, GlobalV::NBANDS);
    std::vector<double> ocp_kb;
    ocp_kb.resize(GlobalV::NBANDS*elecstate->klist->nks-1);
    for (int i = 0; i < ocp_kb.size(); ++i)
    {
        ocp_kb[i] = 1.0;
    }
    testing::internal::CaptureStdout();
    EXPECT_EXIT(elecstate->fixed_weights(ocp_kb), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("size of occupation array is wrong , please check ocp_set"));
}

TEST_F(ElecStateDeathTest,FixedWeightsWarning2)
{
    EXPECT_EQ(GlobalV::NBANDS, 6);
    GlobalV::nelec = 29;
    K_Vectors* klist = new K_Vectors;
    klist->nks = 5;
    elecstate->klist = klist;
    elecstate->wg.create(klist->nks, GlobalV::NBANDS);
    std::vector<double> ocp_kb;
    ocp_kb.resize(GlobalV::NBANDS*elecstate->klist->nks);
    for (int i = 0; i < ocp_kb.size(); ++i)
    {
        ocp_kb[i] = 1.0;
    }
    testing::internal::CaptureStdout();
    EXPECT_EXIT(elecstate->fixed_weights(ocp_kb), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("total number of occupations is wrong , please check ocp_set"));
}

TEST_F(ElecStateTest, CalEBand)
{
    EXPECT_EQ(GlobalV::NBANDS, 6);
    int nks = 5;
    elecstate->wg.create(nks, GlobalV::NBANDS);
    elecstate->ekb.create(nks, GlobalV::NBANDS);
    for (int ik = 0; ik < nks; ++ik)
    {
        for (int ib = 0; ib < GlobalV::NBANDS; ++ib)
        {
            elecstate->ekb(ik, ib) = 1.0;
            elecstate->wg(ik, ib) = 2.0;
        }
    }
    GlobalV::KPAR = 2;
    elecstate->calEBand();
    EXPECT_DOUBLE_EQ(elecstate->f_en.eband, 60.0);
}

TEST_F(ElecStateTest, CalculateWeightsSkipWeights)
{
    EXPECT_FALSE(elecstate->skip_weights);
    elecstate->skip_weights = true;
    EXPECT_NO_THROW(elecstate->calculate_weights());
}

TEST_F(ElecStateDeathTest, CalculateWeightsFixedOccupations)
{
    Occupy::fixed_occupations = true;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(elecstate->calculate_weights(), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("other occupations, not implemented"));
    Occupy::fixed_occupations = false;
}

TEST_F(ElecStateTest, CalculateWeightsIWeights)
{
    EXPECT_FALSE(elecstate->skip_weights);
    int nks = 5;
    K_Vectors* klist = new K_Vectors;
    klist->nks = nks;
    klist->wk.resize(nks);
    for (int ik = 0; ik < nks; ++ik)
    {
        klist->wk[ik] = 2.0;
    }
    klist->isk.resize(nks);
    for (int ik = 0; ik < nks; ++ik)
    {
        klist->isk[ik] = 0;
    }
    elecstate->eferm.ef = 0.0;
    elecstate->klist = klist;
    elecstate->ekb.create(nks, GlobalV::NBANDS);
    for (int ik = 0; ik < nks; ++ik)
    {
        for (int ib = 0; ib < GlobalV::NBANDS; ++ib)
        {
            elecstate->ekb(ik, ib) = 100.0;
        }
    }
    elecstate->wg.create(nks, GlobalV::NBANDS);
    elecstate->calculate_weights();
    EXPECT_DOUBLE_EQ(elecstate->wg(0, 0), 2.0);
    EXPECT_DOUBLE_EQ(elecstate->wg(nks-1, GlobalV::nelec/2-1), 2.0);
    EXPECT_DOUBLE_EQ(elecstate->wg(nks-1, GlobalV::NBANDS-1), 0.0);
    EXPECT_DOUBLE_EQ(elecstate->eferm.ef, 100.0);
    delete klist;
}

TEST_F(ElecStateTest, CalculateWeightsIWeightsTwoFermi)
{
    // get nelec_spin
    GlobalV::TWO_EFERMI = true;
    GlobalV::NSPIN = 2;
    elecstate->init_nelec_spin();
    EXPECT_EQ(elecstate->nelec_spin[0], 5.0);
    EXPECT_EQ(elecstate->nelec_spin[1], 5.0);
    //
    EXPECT_FALSE(elecstate->skip_weights);
    int nks = 5*GlobalV::NSPIN;
    K_Vectors* klist = new K_Vectors;
    klist->nks = nks;
    klist->wk.resize(nks);
    for (int ik = 0; ik < nks; ++ik)
    {
        if(ik<5)
        {
            klist->wk[ik] = 1.1;
        }
        else
        {
            klist->wk[ik] = 1.0;
        }
    }
    klist->isk.resize(nks);
    for (int ik = 0; ik < nks; ++ik)
    {
        if(ik < 5)
        {
            klist->isk[ik] = 0;
        }
        else
        {
            klist->isk[ik] = 1;
        }
    }
    elecstate->eferm.ef_up = 0.0;
    elecstate->eferm.ef_dw = 0.0;
    elecstate->klist = klist;
    elecstate->ekb.create(nks, GlobalV::NBANDS);
    for (int ik = 0; ik < nks; ++ik)
    {
        for (int ib = 0; ib < GlobalV::NBANDS; ++ib)
        {
            if(ik < 5)
            {
                elecstate->ekb(ik, ib) = 100.0;
            }
            else
            {
                elecstate->ekb(ik, ib) = 200.0;
            }
        }
    }
    elecstate->wg.create(nks, GlobalV::NBANDS);
    elecstate->calculate_weights();
    EXPECT_DOUBLE_EQ(elecstate->wg(0, 0), 1.1);
    EXPECT_DOUBLE_EQ(elecstate->wg(nks-1, GlobalV::nelec/2-1), 1.0);
    EXPECT_DOUBLE_EQ(elecstate->wg(nks-1, GlobalV::NBANDS-1), 0.0);
    EXPECT_DOUBLE_EQ(elecstate->eferm.ef_up, 100.0);
    EXPECT_DOUBLE_EQ(elecstate->eferm.ef_dw, 200.0);
    delete klist;
}

TEST_F(ElecStateTest, CalculateWeightsGWeights)
{
    Occupy::use_gaussian_broadening = true;
    EXPECT_FALSE(elecstate->skip_weights);
    int nks = 5;
    K_Vectors* klist = new K_Vectors;
    klist->nks = nks;
    klist->wk.resize(nks);
    for (int ik = 0; ik < nks; ++ik)
    {
        klist->wk[ik] = 2.0;
    }
    klist->isk.resize(nks);
    for (int ik = 0; ik < nks; ++ik)
    {
        klist->isk[ik] = 0;
    }
    elecstate->eferm.ef = 0.0;
    elecstate->klist = klist;
    elecstate->ekb.create(nks, GlobalV::NBANDS);
    for (int ik = 0; ik < nks; ++ik)
    {
        for (int ib = 0; ib < GlobalV::NBANDS; ++ib)
        {
            elecstate->ekb(ik, ib) = 100.0;
        }
    }
    elecstate->wg.create(nks, GlobalV::NBANDS);
    elecstate->calculate_weights();
    // GlobalV::nelec = 10;
    // GlobalV::NBANDS = 6;
    // nks = 5;
    // wg = 10/(5*6) = 0.33333333333
    EXPECT_NEAR(elecstate->wg(0, 0), 0.33333333333, 1e-10);
    EXPECT_NEAR(elecstate->wg(nks-1, GlobalV::nelec/2-1), 0.33333333333, 1e-10);
    EXPECT_NEAR(elecstate->wg(nks-1, GlobalV::NBANDS-1), 0.33333333333,1e-10);
    EXPECT_NEAR(elecstate->eferm.ef, 99.993159296503, 1e-10);
    delete klist;
    Occupy::use_gaussian_broadening = false;
}

TEST_F(ElecStateTest, CalculateWeightsGWeightsTwoFermi)
{
    Occupy::use_gaussian_broadening = true;
    // get nelec_spin
    GlobalV::TWO_EFERMI = true;
    GlobalV::NSPIN = 2;
    elecstate->init_nelec_spin();
    EXPECT_EQ(elecstate->nelec_spin[0], 5.0);
    EXPECT_EQ(elecstate->nelec_spin[1], 5.0);
    //
    EXPECT_FALSE(elecstate->skip_weights);
    int nks = 5*GlobalV::NSPIN;
    K_Vectors* klist = new K_Vectors;
    klist->nks = nks;
    klist->wk.resize(nks);
    for (int ik = 0; ik < nks; ++ik)
    {
        if(ik<5)
        {
            klist->wk[ik] = 1.1;
        }
        else
        {
            klist->wk[ik] = 1.0;
        }
    }
    klist->isk.resize(nks);
    for (int ik = 0; ik < nks; ++ik)
    {
        if(ik < 5)
        {
            klist->isk[ik] = 0;
        }
        else
        {
            klist->isk[ik] = 1;
        }
    }
    elecstate->eferm.ef_up = 0.0;
    elecstate->eferm.ef_dw = 0.0;
    elecstate->klist = klist;
    elecstate->ekb.create(nks, GlobalV::NBANDS);
    for (int ik = 0; ik < nks; ++ik)
    {
        for (int ib = 0; ib < GlobalV::NBANDS; ++ib)
        {
            if(ik < 5)
            {
                elecstate->ekb(ik, ib) = 100.0;
            }
            else
            {
                elecstate->ekb(ik, ib) = 200.0;
            }
        }
    }
    elecstate->wg.create(nks, GlobalV::NBANDS);
    elecstate->calculate_weights();
    // GlobalV::nelec = 10;
    // GlobalV::NBANDS = 6;
    // nks = 10;
    // wg = 10/(10*6) = 0.16666666666
    EXPECT_NEAR(elecstate->wg(0, 0), 0.16666666666, 1e-10);
    EXPECT_NEAR(elecstate->wg(nks-1, GlobalV::nelec/2-1), 0.16666666666, 1e-10);
    EXPECT_NEAR(elecstate->wg(nks-1, GlobalV::NBANDS-1), 0.16666666666, 1e-10);
    EXPECT_NEAR(elecstate->eferm.ef_up, 99.992717105890961, 1e-10);
    EXPECT_NEAR(elecstate->eferm.ef_dw, 199.99315929650351, 1e-10);
    delete klist;
    Occupy::use_gaussian_broadening = false;
}


#undef protected
#undef private