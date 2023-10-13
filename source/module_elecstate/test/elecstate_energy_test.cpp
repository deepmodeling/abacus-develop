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
K_Vectors::K_Vectors()
{
}
K_Vectors::~K_Vectors()
{
}

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
        GlobalV::KPAR = 1;
        GlobalV::NPROC_IN_POOL = 1;
        GlobalV::sc_mag_switch = 1;
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
    // deband_harris + hatree + efiled + gatefield + escon
    EXPECT_DOUBLE_EQ(elecstate->f_en.etot_harris, 0.7);
}

TEST_F(ElecStateEnergyTest, CalEnergiesHarrisImpSol)
{
    elecstate->f_en.deband_harris = 0.1;
    GlobalV::imp_sol = true;
    elecstate->cal_energies(1);
    // deband_harris + hatree + efiled + gatefield + esol_el + esol_cav + escon
    EXPECT_DOUBLE_EQ(elecstate->f_en.etot_harris, 1.6);
}

TEST_F(ElecStateEnergyTest, CalEnergiesHarrisDFTU)
{
    elecstate->f_en.deband_harris = 0.1;
    GlobalV::dft_plus_u = true;
    elecstate->cal_energies(1);
    // deband_harris + hatree + efiled + gatefield + edftu + escon
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
    // deband + hatree + efiled + gatefield + escon
    EXPECT_DOUBLE_EQ(elecstate->f_en.etot, 0.7);
}

TEST_F(ElecStateEnergyTest, CalEnergiesEtotImpSol)
{
    elecstate->f_en.deband = 0.1;
    GlobalV::imp_sol = true;
    elecstate->cal_energies(2);
    // deband + hatree + efiled + gatefield + esol_el + esol_cav + escon
    EXPECT_DOUBLE_EQ(elecstate->f_en.etot, 1.6);
}

TEST_F(ElecStateEnergyTest, CalEnergiesEtotDFTU)
{
    elecstate->f_en.deband = 0.1;
    GlobalV::dft_plus_u = true;
    elecstate->cal_energies(2);
    // deband + hatree + efiled + gatefield + edftu + escon
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

TEST_F(ElecStateEnergyTest, CalBandgapTrivial)
{
    elecstate->cal_bandgap();
    EXPECT_DOUBLE_EQ(elecstate->bandgap, 0.0);
}

TEST_F(ElecStateEnergyTest, CalBandgap)
{
    K_Vectors* klist = new K_Vectors;
    klist->nks = 5;
    elecstate->klist = klist;
    elecstate->ekb.create(klist->nks, GlobalV::NBANDS);
    for (int ik = 0; ik < klist->nks; ik++)
    {
        for (int ib = 0; ib < GlobalV::NBANDS; ib++)
        {
            elecstate->ekb(ik, ib) = ib;
        }
    }
    elecstate->eferm.ef = 2.5;
    elecstate->cal_bandgap();
    EXPECT_DOUBLE_EQ(elecstate->bandgap, 1.0);
}

TEST_F(ElecStateEnergyTest, CalBandgapUpDwTrivial)
{
    elecstate->cal_bandgap_updw();
    EXPECT_DOUBLE_EQ(elecstate->bandgap_up, 0.0);
    EXPECT_DOUBLE_EQ(elecstate->bandgap_dw, 0.0);
}

TEST_F(ElecStateEnergyTest, CalBandgapUpDw)
{
    K_Vectors* klist = new K_Vectors;
    klist->nks = 6;
    klist->isk.resize(6);
    for (int ik = 0; ik < klist->nks; ik++)
    {
        if (ik < 3)
        {
            klist->isk[ik] = 0;
        }
        else
        {
            klist->isk[ik] = 1;
        } 
    }
    elecstate->klist = klist;
    elecstate->ekb.create(klist->nks, GlobalV::NBANDS);
    for (int ik = 0; ik < klist->nks; ik++)
    {
        for (int ib = 0; ib < GlobalV::NBANDS; ib++)
        {
            if (ik < 3)
            {
                elecstate->ekb(ik, ib) = ib;
            }
            else
            {
                elecstate->ekb(ik, ib) = 0.5*ib;
            }
        }
    }
    elecstate->eferm.ef_up = 0.5;
    elecstate->eferm.ef_dw = 2.1;
    elecstate->cal_bandgap_updw();
    EXPECT_DOUBLE_EQ(elecstate->bandgap_up, 1.0);
    EXPECT_DOUBLE_EQ(elecstate->bandgap_dw, 0.5);
}
