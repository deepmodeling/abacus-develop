#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
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
void Charge::init_rho(elecstate::efermi&, ModuleBase::ComplexMatrix const&)
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
 *   - elecstate::ElecState::init_nelec_spin()
 *   - elecstate::ElecState::ElecState()
 *   - elecstate::ElecState::cal_nbands()
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
        GlobalV::TWO_EFERMI = 0.0;
        GlobalV::NBANDS = 6;
        GlobalV::NLOCAL = 6;
        GlobalV::ESOLVER_TYPE = "ksdft";
        GlobalV::LSPINORB = false;
        GlobalV::BASIS_TYPE = "pw";
        GlobalV::md_prec_level = 0;
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
    GlobalV::md_prec_level = 2;
    elecstate->eferm = efermi;
    EXPECT_NO_THROW(elecstate->init_scf(istep, strucfac));
    // delete elecstate->pot is done in the destructor of elecstate
    delete charge;
}
