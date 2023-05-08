#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "module_elecstate/elecstate_getters.h"
#include "module_elecstate/energy.h"

void elecstate::Potential::get_vnew(const Charge* chg, ModuleBase::matrix& vnew)
{
    return;
}

// define mock getter functions
namespace elecstate
{

const int get_rhopw_nrxx()
{
    return 100;
}
const int get_rhopw_nxyz()
{
    return 1000;
}
const double get_ucell_omega()
{
    return 500.0;
}
const double get_ewald_energy()
{
    return 0.1;
}
const double get_hartree_energy()
{
    return 0.2;
}
const double get_etot_efield()
{
    return 0.3;
}
const double get_etot_gatefield()
{
    return 0.4;
}
int xc_functional_type = 3;
const int get_xc_functional_type()
{
    return xc_functional_type;
}
std::string vdw_method = "dftd3";
const std::string get_input_vdw_method()
{
    return vdw_method;
}
#ifdef __LCAO
const double get_dftu_energy()
{
    return 0.5;
}
#endif
#ifdef __DEEPKS
const double get_lcao_deepks_E_delta()
{
    return 0.6;
}
const double get_lcao_deepks_e_delta_band()
{
    return 0.7;
}
#endif
const double get_solvent_model_Ael()
{
    return 0.8;
}
const double get_solvent_model_Acav()
{
    return 0.9;
}
const double get_tot_magnetization()
{
    return 1.1;
}
const double get_abs_magnetization()
{
    return 1.2;
}
const double get_tot_magnetization_nc_x()
{
    return 1.3;
}
const double get_tot_magnetization_nc_y()
{
    return 1.4;
}
const double get_tot_magnetization_nc_z()
{
    return 1.5;
}
#ifdef __EXX
#ifdef __LCAO
const double get_exx_lip_exx_energy()
{
    return 1.6;
}
const bool get_exx_info_ri_real_number()
{
    return true;
}
const double get_exx_lri_double_Eexx()
{
    return 1.7;
}
const std::complex<double> get_exx_lri_complex_Eexx()
{
    return std::complex<double>(1.8, 1.9);
}
const bool get_exx_info_global_cal_exx()
{
    return true;
}
const double get_exx_info_global_hybrid_alpha()
{
    return 2.0;
}
#endif // __LCAO
#endif // __EXX

} // namespace elecstate

/************************************************
 *  unit test of energy.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - energy::energy()
 *   - energy::~energy()
 */

class EnergyTest : public ::testing::Test
{
  protected:
    energy* energy_test;
    void SetUp()
    {
        energy_test = new energy();
    }
    void TearDown()
    {
        delete energy_test;
    }
};

TEST_F(EnergyTest, Energy)
{
    EXPECT_EQ(0, energy_test->etot);
    EXPECT_EQ(0, energy_test->etot_harris);
    EXPECT_EQ(0, energy_test->eband);
    EXPECT_EQ(0, energy_test->deband);
    EXPECT_EQ(0, energy_test->deband_harris);
    EXPECT_EQ(0, energy_test->etxcc);
    EXPECT_EQ(0, energy_test->exx);
    EXPECT_EQ(0, energy_test->demet);
    EXPECT_EQ(0, energy_test->ef);
    EXPECT_EQ(0, energy_test->esol_el);
    EXPECT_EQ(0, energy_test->esol_cav);
}
