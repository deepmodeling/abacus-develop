#include <vector>

#include "gtest/gtest.h"

#define private public
#include "module_elecstate/potentials/potential_new.h"
// mock functions
Structure_Factor::Structure_Factor()
{
}
Structure_Factor::~Structure_Factor()
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
 *   - Constructor: elecstate::Potential(rhopw, rhopw, ucell, vloc, structure_factors, etxc, vtxc) and
 *                  elecstate::Potential::allocate()
 *     - potentials are divided into 2 types: fixed and dynamic
 *     - fixed potentials: loc, gatefield that are independent of rho
 *     - dynamic potentials: hartree, xc, surchem that are dependent of rho
 *   - Getters: elecstate::Potential::get_veff_smooth_data() and elecstate::Potential::get_vofk_smooth_data()
 *     - get the pointers to v_effective and vofk_effective
 *   - PotRegist: elecstate::Potential::pot_egist(components_list)
 *     - add new objects of potentials that are derived classes of PotBase
 *   - CalFixedV: elecstate::Potential::cal_fixed_v()
 *     - calculate the fixed potentials: v_effective_fixed
 *   - CalVeff: elecstate::Potential::cal_v_eff()
 *     - calculate v_effective by adding v_effective_fixed and adding the dynamic potentials
 *   - UpdateFromCharge: elecstate::Potential::update_from_charge()
 *     - calls cal_fixed_v and cal_v_eff to update v_effective from rho
 *   - InitPot: elecstate::Potential::init_pot()
 *     - using istep and update_from_charge to initialize v_effective
 *   - GetVnew: elecstate::Potential::get_vnew()
 *     - used later for scf correction to the forces
 *   - GetEffective: elecstate::Potential::get_effective_v()
 *     - get the matrix reference or double pointer v_effective
 *   - GetEffectiveVOfK: elecstate::Potential::get_effective_vofk()
 *     - get the matrix reference or double pointer vofk_effective
 *   - GetFixedV: elecstate::Potential::get_fixed_v()
 *     - get the double pointer to v_effective_fixed
 */

class PotentialNewTest : public ::testing::Test
{
  protected:
    ModulePW::PW_Basis* rhopw = nullptr;
    ModulePW::PW_Basis* rhodpw = nullptr;
    UnitCell* ucell = nullptr;
    ModuleBase::matrix* vloc = nullptr;
    Structure_Factor* structure_factors = nullptr;
    double* etxc = nullptr;
    double* vtxc = nullptr;
    elecstate::Potential* pot = nullptr;
    virtual void SetUp()
    {
        rhopw = new ModulePW::PW_Basis;
        rhodpw = new ModulePW::PW_Basis;
        ucell = new UnitCell;
        vloc = new ModuleBase::matrix;
        structure_factors = new Structure_Factor();
        etxc = new double;
        vtxc = new double;
        elecstate::Set_GlobalV_Default();
    }
    virtual void TearDown()
    {
        if (rhopw != nullptr)
            delete rhopw;
        if (rhodpw != nullptr)
            delete rhodpw;
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
    pot = new elecstate::Potential(rhopw, rhopw, ucell, vloc, structure_factors, etxc, vtxc);
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
    pot = new elecstate::Potential(rhopw, rhopw, ucell, vloc, structure_factors, etxc, vtxc);
    EXPECT_TRUE(pot->fixed_mode);
    EXPECT_TRUE(pot->dynamic_mode);
    EXPECT_EQ(pot->v_effective_fixed.size(), 100);
    EXPECT_EQ(pot->v_effective.nr, GlobalV::NSPIN);
    EXPECT_EQ(pot->v_effective.nc, 100);
}

TEST_F(PotentialNewTest, ConstructorNRXX0)
{
    rhopw->nrxx = 0;
    pot = new elecstate::Potential(rhopw, rhopw, ucell, vloc, structure_factors, etxc, vtxc);
    EXPECT_TRUE(pot->fixed_mode);
    EXPECT_TRUE(pot->dynamic_mode);
}

TEST_F(PotentialNewTest, ConstructorXC3)
{
    elecstate::tmp_xc_func_type = 3;
    rhopw->nrxx = 100;
    pot = new elecstate::Potential(rhopw, rhopw, ucell, vloc, structure_factors, etxc, vtxc);
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
    pot = new elecstate::Potential(rhopw, rhopw, ucell, vloc, structure_factors, etxc, vtxc);
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
    pot = new elecstate::Potential(rhopw, rhopw, ucell, vloc, structure_factors, etxc, vtxc);
    EXPECT_TRUE(pot->fixed_mode);
    EXPECT_TRUE(pot->dynamic_mode);
    EXPECT_EQ(pot->v_effective_fixed.size(), 100);
    EXPECT_EQ(pot->v_effective.nr, GlobalV::NSPIN);
    EXPECT_EQ(pot->v_effective.nc, 100);
}

TEST_F(PotentialNewTest, Getters)
{
    pot = new elecstate::Potential;
    pot->veff_smooth.create(10, 10);
    pot->vofk_smooth.create(10, 10);
    float* foo;
    foo = pot->get_veff_smooth_data<float>();
    EXPECT_EQ(foo, pot->s_veff_smooth);
    foo = pot->get_vofk_smooth_data<float>();
    EXPECT_EQ(foo, pot->s_vofk_smooth);
    double* doo;
    doo = pot->get_veff_smooth_data<double>();
    EXPECT_EQ(doo, pot->d_veff_smooth);
    doo = pot->get_vofk_smooth_data<double>();
    EXPECT_EQ(doo, pot->d_vofk_smooth);
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
    pot = new elecstate::Potential(rhopw, rhopw, ucell, vloc, structure_factors, etxc, vtxc);
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
    pot = new elecstate::Potential(rhopw, rhopw, ucell, vloc, structure_factors, etxc, vtxc);
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
    pot = new elecstate::Potential(rhopw, rhopw, ucell, vloc, structure_factors, etxc, vtxc);
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
    pot = new elecstate::Potential(rhopw, rhopw, ucell, vloc, structure_factors, etxc, vtxc);
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
    pot = new elecstate::Potential(rhopw, rhopw, ucell, vloc, structure_factors, etxc, vtxc);
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

TEST_F(PotentialNewTest, GetEffectiveVmatrix)
{
    // construct potential
    rhopw->nrxx = 100;
    pot = new elecstate::Potential(rhopw, rhopw, ucell, vloc, structure_factors, etxc, vtxc);
    //
    ModuleBase::matrix v_eff_tmp = pot->get_effective_v();
    const ModuleBase::matrix v_eff_tmp_const = pot->get_effective_v();
    EXPECT_EQ(v_eff_tmp.nr, GlobalV::NSPIN);
    EXPECT_EQ(v_eff_tmp.nc, 100);
    EXPECT_EQ(v_eff_tmp_const.nr, GlobalV::NSPIN);
    EXPECT_EQ(v_eff_tmp_const.nc, 100);
    for (int ir = 0; ir < v_eff_tmp.nr; ir++)
    {
        for (int ic = 0; ic < v_eff_tmp.nc; ic++)
        {
            EXPECT_DOUBLE_EQ(v_eff_tmp(ir, ic), pot->v_effective(ir, ic));
            EXPECT_DOUBLE_EQ(v_eff_tmp_const(ir, ic), pot->v_effective(ir, ic));
        }
    }
}

TEST_F(PotentialNewTest, GetEffectiveVarray)
{
    // construct potential
    rhopw->nrxx = 100;
    pot = new elecstate::Potential(rhopw, rhopw, ucell, vloc, structure_factors, etxc, vtxc);
    //
    double* v_eff_tmp = pot->get_effective_v(0);
    const double* v_eff_tmp_const = pot->get_effective_v(0);
    for (int ic = 0; ic < rhopw->nrxx; ic++)
    {
        EXPECT_DOUBLE_EQ(v_eff_tmp[ic], pot->v_effective(0, ic));
        EXPECT_DOUBLE_EQ(v_eff_tmp_const[ic], pot->v_effective(0, ic));
    }
    v_eff_tmp[0] = 1.0;
    EXPECT_DOUBLE_EQ(pot->v_effective(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(v_eff_tmp_const[0], 1.0);
}

TEST_F(PotentialNewTest, GetEffectiveVarrayNullptr)
{
    pot = new elecstate::Potential;
    EXPECT_EQ(pot->v_effective.nc, 0);
    double* v_eff_tmp = pot->get_effective_v(0);
    const double* v_eff_tmp_const = pot->get_effective_v(0);
    EXPECT_EQ(v_eff_tmp, nullptr);
    EXPECT_EQ(v_eff_tmp_const, nullptr);
}

TEST_F(PotentialNewTest, GetEffectiveVofkmatrix)
{
    // construct potential
    elecstate::tmp_xc_func_type = 3;
    rhopw->nrxx = 100;
    pot = new elecstate::Potential(rhopw, rhopw, ucell, vloc, structure_factors, etxc, vtxc);
    //
    ModuleBase::matrix vofk_eff_tmp = pot->get_effective_vofk();
    const ModuleBase::matrix vofk_eff_tmp_const = pot->get_effective_vofk();
    EXPECT_EQ(vofk_eff_tmp.nr, GlobalV::NSPIN);
    EXPECT_EQ(vofk_eff_tmp.nc, 100);
    EXPECT_EQ(vofk_eff_tmp_const.nr, GlobalV::NSPIN);
    EXPECT_EQ(vofk_eff_tmp_const.nc, 100);
    for (int ir = 0; ir < vofk_eff_tmp.nr; ir++)
    {
        for (int ic = 0; ic < vofk_eff_tmp.nc; ic++)
        {
            EXPECT_DOUBLE_EQ(vofk_eff_tmp(ir, ic), pot->vofk_effective(ir, ic));
            EXPECT_DOUBLE_EQ(vofk_eff_tmp_const(ir, ic), pot->vofk_effective(ir, ic));
        }
    }
}

TEST_F(PotentialNewTest, GetEffectiveVofkarray)
{
    // construct potential
    rhopw->nrxx = 100;
    pot = new elecstate::Potential(rhopw, rhopw, ucell, vloc, structure_factors, etxc, vtxc);
    //
    double* vofk_eff_tmp = pot->get_effective_vofk(0);
    const double* vofk_eff_tmp_const = pot->get_effective_vofk(0);
    for (int ic = 0; ic < rhopw->nrxx; ic++)
    {
        EXPECT_DOUBLE_EQ(vofk_eff_tmp[ic], pot->vofk_effective(0, ic));
        EXPECT_DOUBLE_EQ(vofk_eff_tmp_const[ic], pot->vofk_effective(0, ic));
    }
    vofk_eff_tmp[0] = 1.0;
    EXPECT_DOUBLE_EQ(pot->vofk_effective(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(vofk_eff_tmp_const[0], 1.0);
}

TEST_F(PotentialNewTest, GetEffectiveVofkarrayNullptr)
{
    pot = new elecstate::Potential;
    EXPECT_EQ(pot->v_effective.nc, 0);
    double* vofk_eff_tmp = pot->get_effective_vofk(0);
    const double* vofk_eff_tmp_const = pot->get_effective_vofk(0);
    EXPECT_EQ(vofk_eff_tmp, nullptr);
    EXPECT_EQ(vofk_eff_tmp_const, nullptr);
}

TEST_F(PotentialNewTest, GetFixedV)
{
    rhopw->nrxx = 100;
    pot = new elecstate::Potential(rhopw, rhopw, ucell, vloc, structure_factors, etxc, vtxc);
    EXPECT_TRUE(pot->fixed_mode);
    EXPECT_TRUE(pot->dynamic_mode);
    EXPECT_EQ(pot->v_effective_fixed.size(), 100);
    double* v_eff_fixed_tmp = pot->get_fixed_v();
    const double* v_eff_fixed_tmp_const = pot->get_fixed_v();
    for (int ic = 0; ic < rhopw->nrxx; ic++)
    {
        v_eff_fixed_tmp[ic] = ic;
        EXPECT_DOUBLE_EQ(v_eff_fixed_tmp[ic], pot->v_effective_fixed[ic]);
        EXPECT_DOUBLE_EQ(v_eff_fixed_tmp_const[ic], pot->v_effective_fixed[ic]);
    }
}

TEST_F(PotentialNewTest, GetVeffSmooth)
{
    // construct potential
    rhopw->nrxx = 100;
    elecstate::tmp_xc_func_type = 3;
    pot = new elecstate::Potential(rhopw, rhopw, ucell, vloc, structure_factors, etxc, vtxc);
    //
    ModuleBase::matrix veff_smooth_tmp = pot->get_veff_smooth();
    const ModuleBase::matrix veff_smooth_const_tmp = pot->get_veff_smooth();
    EXPECT_EQ(veff_smooth_tmp.nr, GlobalV::NSPIN);
    EXPECT_EQ(veff_smooth_tmp.nc, 100);
    EXPECT_EQ(veff_smooth_const_tmp.nr, GlobalV::NSPIN);
    EXPECT_EQ(veff_smooth_const_tmp.nc, 100);
    for (int ir = 0; ir < veff_smooth_tmp.nr; ir++)
    {
        for (int ic = 0; ic < veff_smooth_tmp.nc; ic++)
        {
            EXPECT_DOUBLE_EQ(veff_smooth_tmp(ir, ic), pot->veff_smooth(ir, ic));
            EXPECT_DOUBLE_EQ(veff_smooth_const_tmp(ir, ic), pot->veff_smooth(ir, ic));
        }
    }
}

TEST_F(PotentialNewTest, GetVofkSmooth)
{
    // construct potential
    rhopw->nrxx = 100;
    pot = new elecstate::Potential(rhopw, rhopw, ucell, vloc, structure_factors, etxc, vtxc);
    //
    ModuleBase::matrix vofk_smooth_tmp = pot->get_veff_smooth();
    const ModuleBase::matrix vofk_smooth_const_tmp = pot->get_veff_smooth();
    EXPECT_EQ(vofk_smooth_tmp.nr, GlobalV::NSPIN);
    EXPECT_EQ(vofk_smooth_tmp.nc, 100);
    EXPECT_EQ(vofk_smooth_const_tmp.nr, GlobalV::NSPIN);
    EXPECT_EQ(vofk_smooth_const_tmp.nc, 100);
    for (int ir = 0; ir < vofk_smooth_tmp.nr; ir++)
    {
        for (int ic = 0; ic < vofk_smooth_tmp.nc; ic++)
        {
            EXPECT_DOUBLE_EQ(vofk_smooth_tmp(ir, ic), pot->vofk_smooth(ir, ic));
            EXPECT_DOUBLE_EQ(vofk_smooth_const_tmp(ir, ic), pot->vofk_smooth(ir, ic));
        }
    }
}

TEST_F(PotentialNewTest, InterpolateVrsDoubleGrids)
{
    GlobalV::double_grid = true;
    elecstate::tmp_xc_func_type = 3;
    // Init pw_basis
    rhopw->initgrids(4, ModuleBase::Matrix3(1, 0, 0, 0, 1, 0, 0, 0, 1), 4);
    rhopw->initparameters(false, 4);
    rhopw->setuptransform();
    rhopw->collect_local_pw();

    rhodpw->initgrids(4, ModuleBase::Matrix3(1, 0, 0, 0, 1, 0, 0, 0, 1), 6);
    rhodpw->initparameters(false, 6);
    static_cast<ModulePW::PW_Basis_Sup*>(rhodpw)->setuptransform(rhopw);
    rhodpw->collect_local_pw();

    pot = new elecstate::Potential(rhodpw, rhopw, ucell, vloc, structure_factors, etxc, vtxc);

    for (int ir = 0; ir < pot->v_effective.nr; ir++)
    {
        for (int ic = 0; ic < pot->v_effective.nc; ic++)
        {
            pot->v_effective(ir,ic) = ir+ic;
            pot->vofk_effective(ir,ic) = ir+2*ic;
        }
    }

    pot->interpolate_vrs();

    std::vector<double> expect_veff = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
    std::vector<double> expect_vofk =  {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52};

    int index=0;
    for (int ir = 0; ir < pot->veff_smooth.nr; ir++)
    {
        for (int ic = 0; ic < pot->veff_smooth.nc; ic++)
        {
            EXPECT_DOUBLE_EQ(pot->veff_smooth(ir,ic), expect_veff[index]);
            EXPECT_DOUBLE_EQ(pot->vofk_smooth(ir,ic), expect_vofk[index]);
            index++;
        }
    }

}

TEST_F(PotentialNewTest, InterpolateVrsWarningQuit)
{
    GlobalV::double_grid = true;
    // Init pw_basis
    rhopw->initgrids(4, ModuleBase::Matrix3(1, 0, 0, 0, 1, 0, 0, 0, 1), 4);
    rhopw->initparameters(false, 4);
    rhopw->setuptransform();
    rhopw->collect_local_pw();
    rhodpw->gamma_only = false;

    rhodpw->initgrids(4, ModuleBase::Matrix3(1, 0, 0, 0, 1, 0, 0, 0, 1), 6);
    rhodpw->initparameters(false, 6);
    static_cast<ModulePW::PW_Basis_Sup*>(rhodpw)->setuptransform(rhopw);
    rhodpw->collect_local_pw();
    rhodpw->gamma_only = true;

    pot = new elecstate::Potential(rhodpw, rhopw, ucell, vloc, structure_factors, etxc, vtxc);

    EXPECT_EXIT(pot->interpolate_vrs(), ::testing::ExitedWithCode(0), "");
}

TEST_F(PotentialNewTest, InterpolateVrsSingleGrids)
{
    GlobalV::double_grid = false;
    elecstate::tmp_xc_func_type = 3;
    // Init pw_basis
    rhopw->initgrids(4, ModuleBase::Matrix3(1, 0, 0, 0, 1, 0, 0, 0, 1), 4);
    rhopw->initparameters(false, 4);
    rhopw->setuptransform();
    rhopw->collect_local_pw();

    pot = new elecstate::Potential(rhopw, rhopw, ucell, vloc, structure_factors, etxc, vtxc);

    for (int ir = 0; ir < pot->v_effective.nr; ir++)
    {
        for (int ic = 0; ic < pot->v_effective.nc; ic++)
        {
            pot->v_effective(ir,ic) = ir+ic;
            pot->vofk_effective(ir,ic) = ir+2*ic;
        }
    }

    pot->interpolate_vrs();

    for (int ir = 0; ir < pot->veff_smooth.nr; ir++)
    {
        for (int ic = 0; ic < pot->veff_smooth.nc; ic++)
        {
            EXPECT_DOUBLE_EQ(pot->veff_smooth(ir,ic), ir+ic);
            EXPECT_DOUBLE_EQ(pot->vofk_smooth(ir,ic), ir+2*ic);
        }
    }

}

#undef private