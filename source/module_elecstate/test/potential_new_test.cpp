#include <vector>

#include "gtest/gtest.h"

#define private public
#include "module_elecstate/potentials/potential_new.h"
// mock functions
void ModulePW::PW_Basis::initgrids(const double lat0_in,                // unit length (unit in bohr)
                                   const ModuleBase::Matrix3 latvec_in, // Unitcell lattice vectors (unit in lat0)
                                   const double gridecut                // unit in Ry, ecut to set up grids
)
{
}
void ModulePW::PW_Basis::initgrids(const double lat0_in,
                                   const ModuleBase::Matrix3 latvec_in, // Unitcell lattice vectors
                                   const int nx_in,
                                   int ny_in,
                                   int nz_in)
{
}
void ModulePW::PW_Basis::distribute_r()
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
UnitCell::UnitCell()
{
}
UnitCell::~UnitCell()
{
}
Magnetism::Magnetism()
{
}
Magnetism::~Magnetism()
{
}
#ifdef __LCAO
InfoNonlocal::InfoNonlocal()
{
}
InfoNonlocal::~InfoNonlocal()
{
}
#endif
Charge::Charge()
{
}
Charge::~Charge()
{
}

namespace elecstate
{
int tmp_xc_func_type = 1;
int get_xc_func_type()
{
    return tmp_xc_func_type;
}

PotBase* Potential::get_pot_type(const std::string& pot_type)
{
    return new PotBase;
}

void Set_GlobalV_Default()
{
    GlobalV::NSPIN = 1;
    GlobalV::device_flag = "cpu";
    GlobalV::precision_flag = "double";
}
} // namespace elecstate

/************************************************
 *  unit test of potential_new.cpp
 ***********************************************/

/**
 * - Tested Functions:
 */

class PotentialNewTest : public ::testing::Test
{
  protected:
    ModulePW::PW_Basis* rhopw = nullptr;
    UnitCell* ucell = nullptr;
    ModuleBase::matrix* vloc = nullptr;
    ModuleBase::ComplexMatrix* structure_factors = nullptr;
    double* etxc = nullptr;
    double* vtxc = nullptr;
    elecstate::Potential* pot = nullptr;
    virtual void SetUp()
    {
        rhopw = new ModulePW::PW_Basis;
        ucell = new UnitCell;
        vloc = new ModuleBase::matrix;
        structure_factors = new ModuleBase::ComplexMatrix;
        etxc = new double;
        vtxc = new double;
        elecstate::Set_GlobalV_Default();
    }
    virtual void TearDown()
    {
        if (rhopw != nullptr)
            delete rhopw;
        if (ucell != nullptr)
            delete ucell;
        if (vloc != nullptr)
            delete vloc;
        if (structure_factors != nullptr)
            delete structure_factors;
        if (etxc != nullptr)
            delete etxc;
        if (vtxc != nullptr)
            delete vtxc;
        if (pot != nullptr)
            delete pot;
    }
};

TEST_F(PotentialNewTest, ConstructorCPUDouble)
{
    rhopw->nrxx = 100;
    pot = new elecstate::Potential(rhopw, ucell, vloc, structure_factors, etxc, vtxc);
    EXPECT_TRUE(pot->fixed_mode);
    EXPECT_TRUE(pot->dynamic_mode);
    EXPECT_EQ(pot->v_effective_fixed.size(), 100);
    EXPECT_EQ(pot->v_effective.nr, GlobalV::NSPIN);
    EXPECT_EQ(pot->v_effective.nc, 100);
}

TEST_F(PotentialNewTest, ConstructorCPUSingle)
{
    rhopw->nrxx = 100;
    GlobalV::precision_flag = "single";
    pot = new elecstate::Potential(rhopw, ucell, vloc, structure_factors, etxc, vtxc);
    EXPECT_TRUE(pot->fixed_mode);
    EXPECT_TRUE(pot->dynamic_mode);
    EXPECT_EQ(pot->v_effective_fixed.size(), 100);
    EXPECT_EQ(pot->v_effective.nr, GlobalV::NSPIN);
    EXPECT_EQ(pot->v_effective.nc, 100);
}

TEST_F(PotentialNewTest, ConstructorNRXX0)
{
    rhopw->nrxx = 0;
    pot = new elecstate::Potential(rhopw, ucell, vloc, structure_factors, etxc, vtxc);
    EXPECT_TRUE(pot->fixed_mode);
    EXPECT_TRUE(pot->dynamic_mode);
}

TEST_F(PotentialNewTest, ConstructorXC3)
{
    elecstate::tmp_xc_func_type = 3;
    rhopw->nrxx = 100;
    pot = new elecstate::Potential(rhopw, ucell, vloc, structure_factors, etxc, vtxc);
    EXPECT_TRUE(pot->fixed_mode);
    EXPECT_TRUE(pot->dynamic_mode);
    EXPECT_EQ(pot->v_effective_fixed.size(), 100);
    EXPECT_EQ(pot->v_effective.nr, GlobalV::NSPIN);
    EXPECT_EQ(pot->v_effective.nc, 100);
    EXPECT_EQ(pot->vofk_effective.nr, GlobalV::NSPIN);
    EXPECT_EQ(pot->vofk_effective.nc, 100);
}

TEST_F(PotentialNewTest, ConstructorGPUDouble)
{
    // this is just a trivial call to the GPU code
    rhopw->nrxx = 100;
    GlobalV::device_flag = "gpu";
    pot = new elecstate::Potential(rhopw, ucell, vloc, structure_factors, etxc, vtxc);
    EXPECT_TRUE(pot->fixed_mode);
    EXPECT_TRUE(pot->dynamic_mode);
    EXPECT_EQ(pot->v_effective_fixed.size(), 100);
    EXPECT_EQ(pot->v_effective.nr, GlobalV::NSPIN);
    EXPECT_EQ(pot->v_effective.nc, 100);
}

TEST_F(PotentialNewTest, ConstructorGPUSingle)
{
    // this is just a trivial call to the GPU code
    rhopw->nrxx = 100;
    GlobalV::device_flag = "gpu";
    GlobalV::precision_flag = "single";
    pot = new elecstate::Potential(rhopw, ucell, vloc, structure_factors, etxc, vtxc);
    EXPECT_TRUE(pot->fixed_mode);
    EXPECT_TRUE(pot->dynamic_mode);
    EXPECT_EQ(pot->v_effective_fixed.size(), 100);
    EXPECT_EQ(pot->v_effective.nr, GlobalV::NSPIN);
    EXPECT_EQ(pot->v_effective.nc, 100);
}

TEST_F(PotentialNewTest, Getters)
{
    pot = new elecstate::Potential;
    pot->v_effective.create(10, 10);
    pot->vofk_effective.create(10,10);
    float* foo;
    foo = pot->get_v_effective_data<float>();
    EXPECT_EQ(foo, pot->s_v_effective);
    foo = pot->get_vofk_effective_data<float>();
    EXPECT_EQ(foo, pot->s_vofk_effective);
    double* doo;
    doo = pot->get_v_effective_data<double>();
    EXPECT_EQ(doo, pot->d_v_effective);
    doo = pot->get_vofk_effective_data<double>();
    EXPECT_EQ(doo, pot->d_vofk_effective);
    delete foo;
    delete doo;
}

TEST_F(PotentialNewTest, PotRegister)
{
    pot = new elecstate::Potential;
    elecstate::PotBase* pot0 = new elecstate::PotBase;
    pot->components.push_back(pot0);
    EXPECT_EQ(pot->components.size(), 1);
    std::vector<std::string> compnents_list = {"hartree", "xc"};
    pot->pot_register(compnents_list);
    EXPECT_EQ(pot->components.size(), 2);
    EXPECT_FALSE(pot->fixed_done);
}

TEST_F(PotentialNewTest, CalFixedV)
{
    // construct potential
    rhopw->nrxx = 100;
    pot = new elecstate::Potential(rhopw, ucell, vloc, structure_factors, etxc, vtxc);
    EXPECT_TRUE(pot->fixed_mode);
    EXPECT_TRUE(pot->dynamic_mode);
    EXPECT_EQ(pot->v_effective_fixed.size(), 100);
    EXPECT_EQ(pot->v_effective.nr, GlobalV::NSPIN);
    EXPECT_EQ(pot->v_effective.nc, 100);
    //
    std::vector<std::string> compnents_list = {
        "local",
        "hartree",
        "xc",
        "surchem",
        "gatefield"
    };
    std::vector<bool> fixed = {true, false, false, false, true};
    pot->pot_register(compnents_list);
    for (int i = 0; i<compnents_list.size(); i++)
    {
        pot->components[i]->fixed_mode = fixed[i];
    }
    double* vl_pseudo = new double[1000];
    pot->cal_fixed_v(vl_pseudo);
    for (int i = 0; i < pot->v_effective_fixed.size(); i++)
    {
        EXPECT_DOUBLE_EQ(pot->v_effective_fixed[i], 0.0);
    }
    delete[] vl_pseudo;
}

TEST_F(PotentialNewTest, CalVeff)
{
    // construct potential
    rhopw->nrxx = 100;
    pot = new elecstate::Potential(rhopw, ucell, vloc, structure_factors, etxc, vtxc);
    EXPECT_TRUE(pot->fixed_mode);
    EXPECT_TRUE(pot->dynamic_mode);
    EXPECT_EQ(pot->v_effective_fixed.size(), 100);
    EXPECT_EQ(pot->v_effective.nr, GlobalV::NSPIN);
    EXPECT_EQ(pot->v_effective.nc, 100);
    //
    std::vector<std::string> compnents_list = {
        "local",
        "hartree",
        "xc",
        "surchem",
        "gatefield"
    };
    std::vector<bool> dynamic = {false, true, true, true, false};
    pot->pot_register(compnents_list);
    for (int i = 0; i<compnents_list.size(); i++)
    {
        pot->components[i]->dynamic_mode = dynamic[i];
    }
    Charge* chg = new Charge;
    ModuleBase::matrix v_eff;
    v_eff.create(2, 100);
    pot->cal_v_eff(chg,this->ucell,v_eff);
    for (int i = 0; i < pot->v_effective_fixed.size(); i++)
    {
        EXPECT_DOUBLE_EQ(pot->v_effective_fixed[i], 0.0);
    }
    delete chg;
}

TEST_F(PotentialNewTest, UpdateFromCharge)
{
    // construct potential
    rhopw->nrxx = 100;
    pot = new elecstate::Potential(rhopw, ucell, vloc, structure_factors, etxc, vtxc);
    EXPECT_TRUE(pot->fixed_mode);
    EXPECT_TRUE(pot->dynamic_mode);
    EXPECT_EQ(pot->v_effective_fixed.size(), 100);
    EXPECT_EQ(pot->v_effective.nr, GlobalV::NSPIN);
    EXPECT_EQ(pot->v_effective.nc, 100);
    //
    std::vector<std::string> compnents_list = {
        "local",
        "hartree",
        "xc",
        "surchem",
        "gatefield"
    };
    std::vector<bool> fixed = {true, false, false, false, true};
    std::vector<bool> dynamic = {false, true, true, true, false};
    pot->pot_register(compnents_list);
    for (int i = 0; i<compnents_list.size(); i++)
    {
        pot->components[i]->fixed_mode = fixed[i];
        pot->components[i]->dynamic_mode = dynamic[i];
    }
    Charge* chg = new Charge;
    EXPECT_FALSE(pot->fixed_done);
    pot->update_from_charge(chg, this->ucell);
    EXPECT_TRUE(pot->fixed_done);
    delete chg;
}

TEST_F(PotentialNewTest, InitPot)
{
    // construct potential
    rhopw->nrxx = 100;
    pot = new elecstate::Potential(rhopw, ucell, vloc, structure_factors, etxc, vtxc);
    EXPECT_TRUE(pot->fixed_mode);
    EXPECT_TRUE(pot->dynamic_mode);
    EXPECT_EQ(pot->v_effective_fixed.size(), 100);
    EXPECT_EQ(pot->v_effective.nr, GlobalV::NSPIN);
    EXPECT_EQ(pot->v_effective.nc, 100);
    //
    std::vector<std::string> compnents_list = {
        "local",
        "hartree",
        "xc",
        "surchem",
        "gatefield"
    };
    std::vector<bool> fixed = {true, false, false, false, true};
    std::vector<bool> dynamic = {false, true, true, true, false};
    pot->pot_register(compnents_list);
    for (int i = 0; i<compnents_list.size(); i++)
    {
        pot->components[i]->fixed_mode = fixed[i];
        pot->components[i]->dynamic_mode = dynamic[i];
    }
    Charge* chg = new Charge;
    EXPECT_FALSE(pot->fixed_done);
    pot->init_pot(1,chg);
    EXPECT_TRUE(pot->fixed_done);
    delete chg;
}

TEST_F(PotentialNewTest, GetVnew)
{
    // construct potential
    rhopw->nrxx = 100;
    pot = new elecstate::Potential(rhopw, ucell, vloc, structure_factors, etxc, vtxc);
    EXPECT_TRUE(pot->fixed_mode);
    EXPECT_TRUE(pot->dynamic_mode);
    EXPECT_EQ(pot->v_effective_fixed.size(), 100);
    EXPECT_EQ(pot->v_effective.nr, GlobalV::NSPIN);
    EXPECT_EQ(pot->v_effective.nc, 100);
    //
    std::vector<std::string> compnents_list = {
        "local",
        "hartree",
        "xc",
        "surchem",
        "gatefield"
    };
    std::vector<bool> fixed = {true, false, false, false, true};
    std::vector<bool> dynamic = {false, true, true, true, false};
    pot->pot_register(compnents_list);
    for (int i = 0; i<compnents_list.size(); i++)
    {
        pot->components[i]->fixed_mode = fixed[i];
        pot->components[i]->dynamic_mode = dynamic[i];
    }
    Charge* chg = new Charge;
    ModuleBase::matrix vnew;
    pot->get_vnew(chg, vnew);
    EXPECT_EQ(vnew.nr, GlobalV::NSPIN);
    EXPECT_EQ(vnew.nc, 100);
    delete chg;
}

#undef private