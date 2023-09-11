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
int tmp_xc_func_type = 1;
int get_xc_func_type()
{
    return tmp_xc_func_type;
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
void Charge::set_rho_core_paw()
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
    GlobalV::KPAR = 1;
    GlobalV::NPROC_IN_POOL = 1;
}

/************************************************
 *  unit test of elecstate_pw.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - Constructor: elecstate::ElecStatePW constructor and destructor
 *      - including double and single precision versions
 *   - InitRhoData: elecstate::ElecStatePW::init_rho_data()
 *      - get rho and kin_r for ElecStatePW
 *   - ParallelK: elecstate::ElecStatePW::parallelK()
 *      - trivial call due to removing of __MPI
 *   - todo: psiToRho: elecstate::ElecStatePW::psiToRho()
 */

class ElecStatePWTest : public ::testing::Test
{
  protected:
    elecstate::ElecStatePW<double, psi::DEVICE_CPU>* elecstate_pw_d = nullptr;
    elecstate::ElecStatePW<float, psi::DEVICE_CPU>* elecstate_pw_s = nullptr;
    ModulePW::PW_Basis_K* wfcpw = nullptr;
    Charge* chg = nullptr;
    K_Vectors* klist = nullptr;
    ModulePW::PW_Basis* rhopw = nullptr;
    ModulePW::PW_Basis_Big* bigpw = nullptr;
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
        if (elecstate_pw_d != nullptr)
        {
            delete elecstate_pw_d;
        }
        if (elecstate_pw_s != nullptr)
        {
            delete elecstate_pw_s;
        }
    }
};

TEST_F(ElecStatePWTest, ConstructorDouble)
{
    elecstate_pw_d = new elecstate::ElecStatePW<double, psi::DEVICE_CPU>(wfcpw, chg, klist, rhopw, bigpw);
    EXPECT_EQ(elecstate_pw_d->classname, "ElecStatePW");
    EXPECT_EQ(elecstate_pw_d->charge, chg);
    EXPECT_EQ(elecstate_pw_d->klist, klist);
    EXPECT_EQ(elecstate_pw_d->bigpw, bigpw);
}

TEST_F(ElecStatePWTest, ConstructorSingle)
{
    elecstate_pw_s = new elecstate::ElecStatePW<float, psi::DEVICE_CPU>(wfcpw, chg, klist, rhopw, bigpw);
    EXPECT_EQ(elecstate_pw_s->classname, "ElecStatePW");
    EXPECT_EQ(elecstate_pw_s->charge, chg);
    EXPECT_EQ(elecstate_pw_s->klist, klist);
    EXPECT_EQ(elecstate_pw_s->bigpw, bigpw);
}

TEST_F(ElecStatePWTest, InitRhoDataDouble)
{
    elecstate::tmp_xc_func_type = 3;
    chg->nrxx = 1000;
    elecstate_pw_d = new elecstate::ElecStatePW<double, psi::DEVICE_CPU>(wfcpw, chg, klist, rhopw, bigpw);
    elecstate_pw_d->init_rho_data();
    EXPECT_EQ(elecstate_pw_d->init_rho, true);
    EXPECT_EQ(elecstate_pw_d->rho, chg->rho);
    EXPECT_EQ(elecstate_pw_d->kin_r, chg->kin_r);
}

TEST_F(ElecStatePWTest, InitRhoDataSingle)
{
    GlobalV::precision_flag = "single";
    elecstate::tmp_xc_func_type = 3;
    chg->nspin = GlobalV::NSPIN;
    chg->nrxx = 1000;
    elecstate_pw_s = new elecstate::ElecStatePW<float, psi::DEVICE_CPU>(wfcpw, chg, klist, rhopw, bigpw);
    elecstate_pw_s->init_rho_data();
    EXPECT_EQ(elecstate_pw_s->init_rho, true);
    EXPECT_NE(elecstate_pw_s->rho, nullptr);
    EXPECT_NE(elecstate_pw_s->kin_r, nullptr);
}

TEST_F(ElecStatePWTest, ParallelKDouble)
{
    //this is a trivial call due to removing of __MPI
    elecstate_pw_d = new elecstate::ElecStatePW<double, psi::DEVICE_CPU>(wfcpw, chg, klist, rhopw, bigpw);
    EXPECT_NO_THROW(elecstate_pw_d->parallelK());
}

TEST_F(ElecStatePWTest, ParallelKSingle)
{
    //this is a trivial call due to removing of __MPI
    elecstate_pw_s = new elecstate::ElecStatePW<float, psi::DEVICE_CPU>(wfcpw, chg, klist, rhopw, bigpw);
    EXPECT_NO_THROW(elecstate_pw_s->parallelK());
}

#undef protected