#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#define protected public
#include "module_elecstate/elecstate_pw.h"

// mock functions for testing
namespace elecstate
{
double get_ucell_omega()
{
    return 500.0;
}
double get_ucell_tpiba()
{
    return 2.0;
}
int get_xc_func_type()
{
    return 1;
}
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
 *  unit test of elecstate_pw.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - InitNelecSpin:
 */

void Set_GlobalV_Default()
{
    GlobalV::device_flag = "cpu";
    GlobalV::precision_flag = "double";
    GlobalV::DOMAG = false;
    GlobalV::DOMAG_Z = false;
    // Base class dependent
    GlobalV::NSPIN = 1;
    GlobalV::nelec = 10.0;
    GlobalV::nupdown = 0.0;
    GlobalV::TWO_EFERMI = false;
    GlobalV::NBANDS = 6;
    GlobalV::NLOCAL = 6;
    GlobalV::ESOLVER_TYPE = "ksdft";
    GlobalV::LSPINORB = false;
    GlobalV::BASIS_TYPE = "pw";
    GlobalV::md_prec_level = 0;
    GlobalV::KPAR = 1;
    GlobalV::NPROC_IN_POOL = 1;
}

class ElecStatePWTest : public ::testing::Test
{
  protected:
    elecstate::ElecStatePW<double, psi::DEVICE_CPU>* elecstate_pw_d;
    elecstate::ElecStatePW<float, psi::DEVICE_CPU>* elecstate_pw_f;
    ModulePW::PW_Basis_K* wfcpw;
    Charge* chg;
    K_Vectors* klist;
    ModulePW::PW_Basis* rhopw;
    ModulePW::PW_Basis_Big* bigpw;
    void SetUp() override
    {
        Set_GlobalV_Default();
        wfcpw = new ModulePW::PW_Basis_K;
        chg = new Charge;
        klist = new K_Vectors;
        klist->nks = 5;
        rhopw = new ModulePW::PW_Basis;
        bigpw = new ModulePW::PW_Basis_Big;
    }

    void TearDown() override
    {
        delete wfcpw;
        delete chg;
        delete klist;
        delete rhopw;
        delete elecstate_pw_d;
        delete elecstate_pw_f;
    }
};

TEST_F(ElecStatePWTest, Constructor)
{
    // test constructor
    elecstate_pw_d = new elecstate::ElecStatePW<double, psi::DEVICE_CPU>(wfcpw, chg, klist, rhopw, bigpw);
    EXPECT_EQ(elecstate_pw_d->classname, "ElecStatePW");
    EXPECT_EQ(elecstate_pw_d->charge, chg);
    EXPECT_EQ(elecstate_pw_d->klist, klist);
    EXPECT_EQ(elecstate_pw_d->bigpw, bigpw);
    elecstate_pw_f = new elecstate::ElecStatePW<float, psi::DEVICE_CPU>(wfcpw, chg, klist, rhopw, bigpw);
    EXPECT_EQ(elecstate_pw_f->classname, "ElecStatePW");
    EXPECT_EQ(elecstate_pw_f->charge, chg);
    EXPECT_EQ(elecstate_pw_f->klist, klist);
    EXPECT_EQ(elecstate_pw_f->bigpw, bigpw);
}

#undef protected