#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "module_elecstate/elecstate.h"
#include "module_elecstate/elecstate_getters.h"

// mock functions
namespace elecstate
{
int tmp_xc_func_type = 1;
int get_xc_func_type()
{
    return tmp_xc_func_type;
}
void Potential::get_vnew(Charge const*, ModuleBase::matrix&)
{
    return;
}
double ElecState::get_hartree_energy()
{
    return 0.1;
}
double ElecState::get_etot_efield()
{
    return 0.2;
}
double ElecState::get_etot_gatefield()
{
    return 0.3;
}
double ElecState::get_solvent_model_Ael()
{
    return 0.4;
}
double ElecState::get_solvent_model_Acav()
{
    return 0.5;
}
#ifdef __LCAO
double ElecState::get_dftu_energy()
{
    return 0.6;
}
#endif
} // namespace elecstate

/***************************************************************
 *  unit test of functions in elecstate_energy.cpp
 ****************************************************************/

/**
 * - Tested functions:
 */

namespace elecstate
{
class MockElecState : public ElecState
{
  public:
    void Set_GlobalV_Default()
    {
        GlobalV::imp_sol = false;
        GlobalV::dft_plus_u = false;
        // base class
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
};
const double* ElecState::getRho(int spin) const
{
    return &(this->eferm.ef);
} // just for mock
void ElecState::calculate_weights()
{
    return;
} // just for mock
} // namespace elecstate

class ElecStateEnergyTest : public ::testing::Test
{
  protected:
    elecstate::MockElecState* elecstate;
    void SetUp() override
    {
        elecstate = new elecstate::MockElecState;
        elecstate->Set_GlobalV_Default();
    }
    void TearDown() override
    {
        delete elecstate;
    }
};

TEST_F(ElecStateEnergyTest, CalEnergiesHarris)
{
    elecstate->f_en.deband_harris = 0.1;
    elecstate->cal_energies(1);
    // deband_harris + hatree + efiled + gatefield
    EXPECT_DOUBLE_EQ(elecstate->f_en.etot_harris, 0.7);
}

TEST_F(ElecStateEnergyTest, CalEnergiesHarrisImpSol)
{
    elecstate->f_en.deband_harris = 0.1;
    GlobalV::imp_sol = true;
    elecstate->cal_energies(1);
    // deband_harris + hatree + efiled + gatefield + esol_el + esol_cav
    EXPECT_DOUBLE_EQ(elecstate->f_en.etot_harris, 1.6);
}

TEST_F(ElecStateEnergyTest, CalEnergiesHarrisDFTU)
{
    elecstate->f_en.deband_harris = 0.1;
    GlobalV::dft_plus_u = true;
    elecstate->cal_energies(1);
    // deband_harris + hatree + efiled + gatefield + edftu
#ifdef __LCAO
    EXPECT_DOUBLE_EQ(elecstate->f_en.etot_harris, 1.3);
#else
    EXPECT_DOUBLE_EQ(elecstate->f_en.etot_harris, 0.7);
#endif
}

TEST_F(ElecStateEnergyTest, CalEnergiesEtot)
{
    elecstate->f_en.deband = 0.1;
    elecstate->cal_energies(2);
    // deband + hatree + efiled + gatefield
    EXPECT_DOUBLE_EQ(elecstate->f_en.etot, 0.7);
}

TEST_F(ElecStateEnergyTest, CalEnergiesEtotImpSol)
{
    elecstate->f_en.deband = 0.1;
    GlobalV::imp_sol = true;
    elecstate->cal_energies(2);
    // deband + hatree + efiled + gatefield + esol_el + esol_cav
    EXPECT_DOUBLE_EQ(elecstate->f_en.etot, 1.6);
}

TEST_F(ElecStateEnergyTest, CalEnergiesEtotDFTU)
{
    elecstate->f_en.deband = 0.1;
    GlobalV::dft_plus_u = true;
    elecstate->cal_energies(2);
    // deband + hatree + efiled + gatefield + edftu
#ifdef __LCAO
    EXPECT_DOUBLE_EQ(elecstate->f_en.etot, 1.3);
#else
    EXPECT_DOUBLE_EQ(elecstate->f_en.etot, 0.7);
#endif
}

TEST_F(ElecStateEnergyTest, CalConverged)
{
    elecstate->cal_converged();
    EXPECT_TRUE(elecstate->vnew_exist);
    EXPECT_DOUBLE_EQ(elecstate->f_en.descf, 0.0);
}