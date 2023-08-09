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
 *   - Constructor: elecstate::Potential(rhopw, ucell, vloc, structure_factors, etxc, vtxc) and elecstate::Potential::allocate()
 *     - potentials are divided into 2 types: fixed and dynamic
 *     - fixed potentials: loc, gatefield that are independent of rho
 *     - dynamic potentials: hartree, xc, surchem that are dependent of rho
 *   - Getters: elecstate::Potential::get_v_effective_data() and elecstate::Potential::get_vofk_effective_data()
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
    UnitCell* ucell = nullptr;
    ModuleBase::matrix* vloc = nullptr;
    Structure_Factor* structure_factors = nullptr;
    double* etxc = nullptr;
    double* vtxc = nullptr;
    elecstate::Potential* pot = nullptr;
    virtual void SetUp()
    {
        rhopw = new ModulePW::PW_Basis;
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
    pot = new elecstate::Potential(rhopw, ucell, vloc, structure_factors, etxc, vtxc);
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
    pot = new elecstate::Potential(rhopw, ucell, vloc, structure_factors, etxc, vtxc);
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
    pot = new elecstate::Potential(rhopw, ucell, vloc, structure_factors, etxc, vtxc);
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
    pot = new elecstate::Potential(rhopw, ucell, vloc, structure_factors, etxc, vtxc);
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
    pot = new elecstate::Potential(rhopw, ucell, vloc, structure_factors, etxc, vtxc);
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
#undef private

#include "mpi.h"
int main(int argc, char **argv)
{
        MPI_Init(&argc, &argv);
        testing::InitGoogleTest(&argc, argv);

        int nproc;
        int my_rank;

        MPI_Comm_size(MPI_COMM_WORLD,&nproc);
        MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

        int result = RUN_ALL_TESTS();

        MPI_Finalize();

        return result;
}